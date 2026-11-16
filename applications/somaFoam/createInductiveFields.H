// Inductive (ICP) option. electroMagnetics: dischargeType CCP | ICP (default CCP).
// For ICP, the azimuthal vector potential A = AR + i AI is solved on a second mesh
// region "em" (plasma cells first, then window and coil space), with the real and
// imaginary parts block-coupled through the plasma conductivity. The cycle-averaged
// inductive power density Pind is registered on the plasma mesh and added to the
// electron energy equation by efullImplicit.

const word dischargeType(electroMagnetics.lookupOrDefault<word>("dischargeType", "CCP"));
if (dischargeType != "CCP" && dischargeType != "ICP")
{
    FatalError << "in electroMagnetics, unknown dischargeType " << dischargeType
        << nl << "valid: CCP ICP" << abort(FatalError);
}
const bool isICP = (dischargeType == "ICP");

autoPtr<fvMesh> emMeshPtr;
autoPtr<volScalarField> ARPtr, AIPtr, PindPtr;
autoPtr<volVector2Field> APtr;
scalarField JcoilEM;
scalarField icpPindPrev;
labelList emCellOfPlasma;
scalar icpFreq = 0, icpPower = 0, icpWedge = 360, icpCoilI = 0, icpRamp = 0, icpRelax = 1;
label icpInterval = 1;
// adaptive update: the field is re-solved when the conductivity change since the last
// solve, weighted by the field intensity, exceeds updateTolerance, or after
// maxUpdateInterval steps; the deposited power is rescaled to the set power in between
scalar icpTol = 0, icpPsetLast = 0, icpP1 = 0;
label icpMaxInterval = 1, icpLastSolve = -1, icpNSolves = 0;
scalarField icpSigRLast, icpSigILast;
bool icpPending = true;
bool icpVerify = false;
scalar icpVn0 = 0, icpVnu = 0, icpVR = 0, icpVL = 0;
label icpEIndex = -1;

if (isICP)
{
    const dictionary& icpDict = electroMagnetics.subDict("ICPCoeffs");
    icpFreq = readScalar(icpDict.lookup("frequency"));
    icpPower = readScalar(icpDict.lookup("power"));
    icpWedge = readScalar(icpDict.lookup("wedgeAngle"));
    icpInterval = readLabel(icpDict.lookup("updateInterval"));
    icpTol = icpDict.lookupOrDefault<scalar>("updateTolerance", 0.0);
    icpMaxInterval = icpDict.lookupOrDefault<label>("maxUpdateInterval", icpInterval);
    // power ramp time (s) and under-relaxation of the deposited power, as in the
    // iterative module coupling of hybrid plasma equipment models
    icpRamp = icpDict.lookupOrDefault<scalar>("powerRampTime", 0.0);
    icpRelax = icpDict.lookupOrDefault<scalar>("powerRelaxation", 1.0);
    const scalarList coilR(icpDict.lookup("coilR"));
    const scalar coilZ = readScalar(icpDict.lookup("coilZ"));
    const scalar coilW = readScalar(icpDict.lookup("coilW"));
    icpEIndex = composition.species()["electron"];

    emMeshPtr.reset
    (
        new fvMesh(IOobject("em", runTime.timeName(), runTime, IOobject::MUST_READ))
    );
    const fvMesh& emMesh = emMeshPtr();

    ARPtr.reset(new volScalarField(IOobject("AR", runTime.timeName(), emMesh,
        IOobject::MUST_READ, IOobject::AUTO_WRITE), emMesh));
    AIPtr.reset(new volScalarField(IOobject("AI", runTime.timeName(), emMesh,
        IOobject::MUST_READ, IOobject::AUTO_WRITE), emMesh));
    // Block field (AR, AI) built from the boundary types of AR, so that it needs no
    // file and is created on every processor without decomposing a vector2 field
    wordList aTypes(ARPtr().boundaryField().size());
    forAll(aTypes, p)
    {
        const fvPatchScalarField& pf = ARPtr().boundaryField()[p];
        const word& t = pf.type();
        aTypes[p] = (t == "wedge" || t == "empty" || t == "processor" || t == "symmetryPlane")
            ? t : (pf.fixesValue() ? word("fixedValue") : word("zeroGradient"));
    }
    APtr.reset(new volVector2Field(IOobject("A", runTime.timeName(), emMesh,
        IOobject::NO_READ, IOobject::NO_WRITE), emMesh,
        dimensioned<vector2>("zero", ARPtr().dimensions(), vector2::zero), aTypes));

    PindPtr.reset(new volScalarField(IOobject("Pind", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::AUTO_WRITE), mesh,
        dimensionedScalar("zero", dimensionSet(1, -1, -3, 0, 0), 0.0)));

    // coil current density for 1 A per turn: each square conductor of side coilW is
    // distributed over the cells by the exact overlap of its (r, z) cross-section with
    // the (r, z) extent of each cell, so every turn carries 1 A on any mesh
    const vectorField& Cem = emMesh.C().internalField();
    const pointField& pts = emMesh.points();
    const labelListList& cPts = emMesh.cellPoints();
    JcoilEM.setSize(Cem.size(), 0.0);
    scalar coilCurrentCheck = 0;
    forAll(Cem, i)
    {
        scalar r0 = GREAT, r1 = -GREAT, z0 = GREAT, z1 = -GREAT;
        forAll(cPts[i], k)
        {
            const point& p = pts[cPts[i][k]];
            const scalar r = Foam::sqrt(sqr(p.x()) + sqr(p.y()));
            r0 = min(r0, r); r1 = max(r1, r); z0 = min(z0, p.z()); z1 = max(z1, p.z());
        }
        forAll(coilR, c)
        {
            const scalar lr = max(0.0, min(r1, coilR[c] + coilW/2) - max(r0, coilR[c] - coilW/2));
            const scalar lz = max(0.0, min(z1, coilZ + coilW/2) - max(z0, coilZ - coilW/2));
            if (lr*lz > 0)
            {
                JcoilEM[i] += lr*lz/sqr(coilW)/((r1 - r0)*(z1 - z0));
                coilCurrentCheck += lr*lz/sqr(coilW);
            }
        }
    }
    coilCurrentCheck = returnReduce(coilCurrentCheck, sumOp<scalar>());

    if (icpDict.found("verification"))
    {
        const dictionary& vd = icpDict.subDict("verification");
        icpVerify = true;
        icpVn0 = readScalar(vd.lookup("n0"));
        icpVnu = readScalar(vd.lookup("nuM"));
        icpVR = readScalar(vd.lookup("R"));
        icpVL = readScalar(vd.lookup("L"));
    }

    Info<< "ICP: em mesh " << returnReduce(Cem.size(), sumOp<label>())
        << " cells, coil cells " << returnReduce(label(sum(pos(JcoilEM - SMALL))), sumOp<label>())
        << ", total coil current " << coilCurrentCheck << " A for " << coilR.size()
        << " turns (expected " << coilR.size() << " A)" << endl;
}
