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
// verification profile: n0 J0(2.405 r/R) sin(pi z/L) (besselSine) or n0 for
// r < R and z < L (uniform)
bool icpVUniform = false;
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
    // turns of the coil: coilTurns ((r z) ...), or coilR (r ...) at coilZ
    const List<vector2> coilTurns(inductiveField::coilTurns(icpDict));
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

    // coil current density for 1 A per turn, by the exact overlap of each
    // square conductor with the cells (inductiveField.H)
    const vectorField& Cem = emMesh.C().internalField();
    const scalar coilCurrentCheck =
        inductiveField::coilCurrentDensity(emMesh, coilTurns, coilW, JcoilEM);

    if (icpDict.found("verification"))
    {
        const dictionary& vd = icpDict.subDict("verification");
        icpVerify = true;
        icpVn0 = readScalar(vd.lookup("n0"));
        icpVnu = readScalar(vd.lookup("nuM"));
        icpVR = readScalar(vd.lookup("R"));
        icpVL = readScalar(vd.lookup("L"));
        icpVUniform =
            (vd.lookupOrDefault<word>("profile", "besselSine") == "uniform");
    }

    Info<< "ICP: em mesh " << returnReduce(Cem.size(), sumOp<label>())
        << " cells, coil cells " << returnReduce(label(sum(pos(JcoilEM - SMALL))), sumOp<label>())
        << ", total coil current " << coilCurrentCheck << " A for " << coilTurns.size()
        << " turns (expected " << coilTurns.size() << " A)" << endl;
}
