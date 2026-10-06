/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

Application
    ionTracker

Description
    Ion energy and angular distributions at the walls from the fields of a
    fluid solution: a test-particle Monte Carlo calculation run after
    somaFoam (in the manner of the fluid/Monte Carlo and PCMCM approaches
    of Economou and of Kushner and co-workers).

    Ions are launched where the fluid solution produces them, at random
    phases of the period, with the velocity of a gas atom. They are moved
    in the electric field of the fluid solution, E(x, t), which is read at
    nPhases times of one period and interpolated linearly in space (from
    the cell value and gradient) and in time, and they collide with the
    background gas (null-collision method): charge exchange, in which the
    ion takes the velocity of the atom, and isotropic scattering in the
    centre-of-mass frame. An ion that reaches a patch is scored with its
    energy, its angle to the patch normal and the phase of arrival.

    system/ionTrackerDict:

    \verbatim
    ion             Arp1;       // flux field J_<ion> gives the source
    ionMass         39.948;     // [amu]; gasMass likewise, default ionMass
    charge          1;

    electricField   E;
    gasDensity      Ar;         // number density of the background gas
    gasTemperature  T;

    startTime       2.975e-6;   // one period of fields written by somaFoam
    period          2.5e-8;
    nPhases         50;

    source          fluxDivergence;   // period mean of div(J_<ion>), where
                                      // positive; or density (ion density)
                                      // or the name of a field
    nParticles      20000;
    patches         (electrode ground);

    crossSections   PhelpsAr;   // Ar+ in Ar (Phelps 1994), or
    // crossSections constant;  chargeExchange 5e-19;  isotropic 0;  [m2]
    // crossSections LXCat;     // tables from a file in the LXCat format
    //     file "constant/Ar+_Ar.txt";
    //     energyFrame centreOfMass;   // or laboratory: the frame of the
    //                                 // energies of the file (required)
    //     backwardProcess BACKSCAT;   // keyword of the block used for
    //                                 // charge exchange (default)
    //     isotropicProcess ISOTROPIC; // isotropic scattering (default);
    //                                 // none: no such block
    //     species "Ar+";              // optional: text that the line
    //                                 // after the keyword must contain
    // crossSections table;     // tables of (energy [eV], cross section [m2])
    //     energyFrame laboratory;
    //     chargeExchangeTable ((0.01 1e-18) (100 5e-19));
    //     isotropicTable ((0.01 2e-18) (100 1e-20));   // optional
    maxCollisionEnergy 200;     // [eV] for the null-collision frequency

    maxTime         1e-3;       // longest track [s]
    stepsPerPeriod  200;        // at least this many steps per period
    cellFraction    0.3;        // largest step as a fraction of the cell
    nEnergyBins     200;
    nAngleBins      90;
    // maxEnergy    50;         // [eV]; default: the largest energy scored
    writeParticles  no;         // list of all scored ions
    seed            1234;
    \endverbatim

    Output in <case>/ionTracker: for each patch the energy distribution
    (normalised to unit integral), the angular distribution and the joint
    distribution, and a summary with the fluxes and mean energies. The time
    directory of startTime receives ionTrackerDensity, ionTrackerEnergy and
    ionTrackerVelocity: the ion density, mean energy [eV] and mean velocity
    of the tracked ions, for comparison with the fluid solution.

    Serial runs on a static mesh; all positive ions are treated as the one
    species given.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "ionParticle.H"
#include "Random.H"
#include "OFstream.H"
#include "IFstream.H"
#include "IStringStream.H"
#include "Tuple2.H"
#include "mathematicalConstants.H"

namespace Foam
{
    defineParticleTypeNameAndDebug(ionParticle, 0);
    defineTemplateTypeNameAndDebug(Cloud<ionParticle>, 0);
}

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Cross sections [m2] for charge exchange and isotropic scattering as
//  functions of the energy eps [eV] of the ion relative to the atom at rest
class crossSectionModel
{
    word type_;
    scalar chargeExchange_;
    scalar isotropic_;

    //- Tabulated cross sections: ln(energy) and ln(cross section)
    scalarField isoLogE_;
    scalarField isoLogS_;
    scalarField cxLogE_;
    scalarField cxLogS_;

    //- Factor from the energy of the ion relative to an atom at rest to
    //  the energy of the tables (1, or the centre-of-mass fraction)
    scalar energyFactor_;


    //- Store a table of (energy, cross section) pairs as logarithms
    static void setTable
    (
        const List<Tuple2<scalar, scalar> >& table,
        scalarField& logE,
        scalarField& logS
    )
    {
        logE.setSize(table.size());
        logS.setSize(table.size());

        forAll(table, i)
        {
            logE[i] = Foam::log(max(table[i].first(), 1e-30));
            logS[i] = Foam::log(max(table[i].second(), 1e-60));

            if (i > 0 && logE[i] <= logE[i - 1])
            {
                FatalErrorIn("crossSectionModel::setTable")
                    << "Energies of a cross section table must increase"
                    << exit(FatalError);
            }
        }
    }

    //- Interpolate a table linearly in the logarithms; constant beyond
    //  its ends
    static scalar interpolate
    (
        const scalarField& logE,
        const scalarField& logS,
        const scalar eps
    )
    {
        if (logE.empty())
        {
            return 0;
        }

        const scalar x = Foam::log(max(eps, 1e-30));
        const label n = logE.size();

        if (x <= logE[0])
        {
            return Foam::exp(logS[0]);
        }

        if (x >= logE[n - 1])
        {
            return Foam::exp(logS[n - 1]);
        }

        label low = 0;
        label high = n - 1;

        while (high - low > 1)
        {
            const label mid = (low + high)/2;

            if (logE[mid] <= x)
            {
                low = mid;
            }
            else
            {
                high = mid;
            }
        }

        const scalar w = (x - logE[low])/(logE[high] - logE[low]);

        return Foam::exp((1 - w)*logS[low] + w*logS[high]);
    }

    //- Read the block of a file in the LXCat format that starts with a
    //  line equal to the keyword and, if a species is given, has it in its
    //  next line: the table between the two lines of dashes
    static List<Tuple2<scalar, scalar> > readLXCat
    (
        const fileName& file,
        const word& keyword,
        const string& species
    )
    {
        IFstream is(file);

        if (!is.good())
        {
            FatalErrorIn("crossSectionModel::readLXCat")
                << "Cannot read " << file << exit(FatalError);
        }

        DynamicList<Tuple2<scalar, scalar> > table;

        // 0: looking for the keyword, 1: line after the keyword,
        // 2: looking for the first dashes, 3: in the table
        label state = 0;

        while (is.good())
        {
            string line;
            is.getLine(line);

            // Trim
            const string::size_type first = line.find_first_not_of(" \t\r");
            const string::size_type last = line.find_last_not_of(" \t\r");

            const string trimmed
            (
                first == string::npos
              ? string("")
              : string(line.substr(first, last - first + 1))
            );

            if (state == 0)
            {
                if (trimmed == keyword)
                {
                    state = 1;
                }
            }
            else if (state == 1)
            {
                if (species.empty() || trimmed.find(species) != string::npos)
                {
                    state = 2;
                }
                else
                {
                    state = 0;
                }
            }
            else if (state == 2)
            {
                if (trimmed.size() >= 5 && trimmed.substr(0, 5) == "-----")
                {
                    state = 3;
                }
            }
            else
            {
                if (trimmed.size() >= 5 && trimmed.substr(0, 5) == "-----")
                {
                    break;
                }

                IStringStream row(trimmed);
                scalar e, sig;
                row >> e >> sig;

                table.append(Tuple2<scalar, scalar>(e, sig));
            }
        }

        if (table.size() < 2)
        {
            FatalErrorIn("crossSectionModel::readLXCat")
                << "No table for the process " << keyword
                << (species.empty() ? string("") : string(" with ") + species)
                << " in " << file << exit(FatalError);
        }

        Info<< "    " << keyword << ": " << table.size()
            << " points from " << table[0].first() << " to "
            << table[table.size() - 1].first() << " eV" << endl;

        List<Tuple2<scalar, scalar> > result;
        result.transfer(table);

        return result;
    }

public:

    crossSectionModel
    (
        const dictionary& dict,
        const fileName& casePath,
        const scalar ionMass,
        const scalar gasMass
    )
    :
        type_(dict.lookup("crossSections")),
        chargeExchange_(0),
        isotropic_(0),
        energyFactor_(1)
    {
        if (type_ == "constant")
        {
            chargeExchange_ = readScalar(dict.lookup("chargeExchange"));
            isotropic_ = dict.lookupOrDefault<scalar>("isotropic", 0);
        }
        else if (type_ == "LXCat" || type_ == "table")
        {
            // The frame of the energies of the tables must be stated
            const word frame(dict.lookup("energyFrame"));

            if (frame == "centreOfMass")
            {
                energyFactor_ = gasMass/(ionMass + gasMass);
            }
            else if (frame != "laboratory")
            {
                FatalIOErrorIn("crossSectionModel", dict)
                    << "energyFrame must be laboratory (energy of the ion"
                    << " relative to an atom at rest) or centreOfMass"
                    << exit(FatalIOError);
            }

            if (type_ == "LXCat")
            {
                fileName file(dict.lookup("file"));
                file.expand();

                if (file[0] != '/')
                {
                    file = casePath/file;
                }

                const string species
                (
                    dict.lookupOrDefault<string>("species", "")
                );

                Info<< "Cross sections from " << file << endl;

                setTable
                (
                    readLXCat
                    (
                        file,
                        dict.lookupOrDefault<word>
                        (
                            "backwardProcess",
                            "BACKSCAT"
                        ),
                        species
                    ),
                    cxLogE_,
                    cxLogS_
                );

                const word isoKeyword
                (
                    dict.lookupOrDefault<word>("isotropicProcess", "ISOTROPIC")
                );

                if (isoKeyword != "none")
                {
                    setTable
                    (
                        readLXCat(file, isoKeyword, species),
                        isoLogE_,
                        isoLogS_
                    );
                }
            }
            else
            {
                setTable
                (
                    List<Tuple2<scalar, scalar> >
                    (
                        dict.lookup("chargeExchangeTable")
                    ),
                    cxLogE_,
                    cxLogS_
                );

                if (dict.found("isotropicTable"))
                {
                    setTable
                    (
                        List<Tuple2<scalar, scalar> >
                        (
                            dict.lookup("isotropicTable")
                        ),
                        isoLogE_,
                        isoLogS_
                    );
                }
            }
        }
        else if (type_ != "PhelpsAr")
        {
            FatalIOErrorIn("crossSectionModel", dict)
                << "Unknown crossSections " << type_
                << "; valid entries are PhelpsAr, constant, LXCat and table"
                << exit(FatalIOError);
        }
    }

    //- Isotropic part
    scalar isotropic(const scalar eps) const
    {
        if (type_ == "constant")
        {
            return isotropic_;
        }
        else if (type_ != "PhelpsAr")
        {
            return interpolate(isoLogE_, isoLogS_, energyFactor_*eps);
        }

        // Phelps, J. Appl. Phys. 76 (1994) 747, Ar+ in Ar
        const scalar e = max(eps, 1e-4);

        return
            2e-19/(Foam::sqrt(e)*(1 + e))
          + 3e-19*e/Foam::pow(1 + e/3, 2.3);
    }

    //- Backward part (charge exchange)
    scalar chargeExchange(const scalar eps) const
    {
        if (type_ == "constant")
        {
            return chargeExchange_;
        }
        else if (type_ != "PhelpsAr")
        {
            return interpolate(cxLogE_, cxLogS_, energyFactor_*eps);
        }

        const scalar e = max(eps, 1e-4);

        const scalar momentum =
            1.15e-18*Foam::pow(e, -0.1)*Foam::pow(1 + 0.015/e, 0.6);

        return max(0.5*(momentum - isotropic(e)), 0.0);
    }
};


//- Random unit vector
vector randomDirection(Random& rndGen)
{
    const scalar cosTheta = 2*rndGen.scalar01() - 1;
    const scalar sinTheta = Foam::sqrt(max(1 - sqr(cosTheta), 0.0));
    const scalar phi = 2*mathematicalConstant::pi*rndGen.scalar01();

    return vector(sinTheta*Foam::cos(phi), sinTheta*Foam::sin(phi), cosTheta);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    if (Pstream::parRun())
    {
        FatalErrorIn(args.executable())
            << "ionTracker runs in serial" << exit(FatalError);
    }

    const scalar eCharge = 1.602176634e-19;
    const scalar amu = 1.66053906660e-27;
    const scalar kB = 1.380649e-23;

    IOdictionary dict
    (
        IOobject
        (
            "ionTrackerDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const word ionName(dict.lookup("ion"));
    const scalar ionMass = readScalar(dict.lookup("ionMass"))*amu;
    const scalar gasMass =
        dict.lookupOrDefault<scalar>("gasMass", ionMass/amu)*amu;
    const scalar charge = dict.lookupOrDefault<scalar>("charge", 1)*eCharge;

    const word EName(dict.lookupOrDefault<word>("electricField", "E"));
    const word gasDensityName(dict.lookup("gasDensity"));
    const word gasTemperatureName
    (
        dict.lookupOrDefault<word>("gasTemperature", "T")
    );

    const scalar startTime = readScalar(dict.lookup("startTime"));
    const scalar period = readScalar(dict.lookup("period"));
    const label nPhases = readLabel(dict.lookup("nPhases"));

    const word sourceType
    (
        dict.lookupOrDefault<word>("source", "fluxDivergence")
    );
    const label nParticles = readLabel(dict.lookup("nParticles"));
    const wordList patchNames(dict.lookup("patches"));

    const crossSectionModel sigma(dict, runTime.path(), ionMass, gasMass);
    const scalar maxCollisionEnergy =
        dict.lookupOrDefault<scalar>("maxCollisionEnergy", 200);

    const scalar maxTime = dict.lookupOrDefault<scalar>("maxTime", 1e-3);
    const label stepsPerPeriod =
        dict.lookupOrDefault<label>("stepsPerPeriod", 200);
    const scalar cellFraction =
        dict.lookupOrDefault<scalar>("cellFraction", 0.3);
    const label nEnergyBins = dict.lookupOrDefault<label>("nEnergyBins", 200);
    const label nAngleBins = dict.lookupOrDefault<label>("nAngleBins", 90);
    scalar maxEnergy = dict.lookupOrDefault<scalar>("maxEnergy", -1);
    const Switch writeParticles
    (
        dict.lookupOrDefault<Switch>("writeParticles", false)
    );

    Random rndGen(dict.lookupOrDefault<label>("seed", 1234));


    // Fields of one period

    const instantList times = runTime.times();

    PtrList<volVectorField> EFields(nPhases);
    PtrList<volTensorField> gradEFields(nPhases);

    scalarField source(mesh.nCells(), 0.0);

    Info<< "Reading " << EName << " at " << nPhases << " phases of the period"
        << " from " << startTime << " s" << endl;

    forAll(EFields, phaseI)
    {
        const scalar t = startTime + phaseI*period/nPhases;

        label best = -1;
        scalar bestDiff = GREAT;

        forAll(times, i)
        {
            if (times[i].name() != "constant")
            {
                const scalar diff = mag(times[i].value() - t);

                if (diff < bestDiff)
                {
                    bestDiff = diff;
                    best = i;
                }
            }
        }

        if (best < 0 || bestDiff > 0.25*period/nPhases)
        {
            FatalErrorIn(args.executable())
                << "No time directory near " << t << " s (phase " << phaseI
                << "). Write the fields " << nPhases << " times per period"
                << " over one period, from startTime." << exit(FatalError);
        }

        runTime.setTime(times[best], best);

        EFields.set
        (
            phaseI,
            new volVectorField
            (
                IOobject
                (
                    EName,
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh
            )
        );

        gradEFields.set
        (
            phaseI,
            new volTensorField
            (
                IOobject
                (
                    "grad" + EName,
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                fvc::grad(EFields[phaseI])
            )
        );

        if (sourceType == "fluxDivergence")
        {
            const volVectorField J
            (
                IOobject
                (
                    "J_" + ionName,
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh
            );

            source += fvc::div(J)().internalField()/nPhases;
        }
        else
        {
            const volScalarField s
            (
                IOobject
                (
                    sourceType == "density" ? ionName : sourceType,
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh
            );

            source += s.internalField()/nPhases;
        }
    }

    // Back to the first phase for the gas fields and for the output
    {
        label first = 0;
        scalar bestDiff = GREAT;

        forAll(times, i)
        {
            if
            (
                times[i].name() != "constant"
             && mag(times[i].value() - startTime) < bestDiff
            )
            {
                bestDiff = mag(times[i].value() - startTime);
                first = i;
            }
        }

        runTime.setTime(times[first], first);
    }

    const volScalarField gasDensity
    (
        IOobject
        (
            gasDensityName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh
    );

    const volScalarField gasTemperature
    (
        IOobject
        (
            gasTemperatureName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh
    );


    // Source: ions per second in each cell, and its cumulative distribution

    const scalarField& V = mesh.V();

    source = max(source, scalar(0));

    scalarField cumulative(mesh.nCells());
    scalar ionsPerSecond = 0;

    forAll(source, cellI)
    {
        ionsPerSecond += source[cellI]*V[cellI];
        cumulative[cellI] = ionsPerSecond;
    }

    if (ionsPerSecond <= 0)
    {
        FatalErrorIn(args.executable())
            << "The source " << sourceType << " is nowhere positive"
            << exit(FatalError);
    }

    cumulative /= ionsPerSecond;

    // What one tracked ion stands for [1/s] (for source fluxDivergence)
    const scalar weight = ionsPerSecond/nParticles;


    // Geometry

    const Vector<label>& directions = mesh.geometricD();

    // Smallest dimension of each cell
    scalarField cellSize(mesh.nCells());
    {
        const cellList& cells = mesh.cells();
        const scalarField& magSf = mesh.magSf().internalField();

        forAll(cells, cellI)
        {
            scalar maxArea = 0;

            forAll(cells[cellI], i)
            {
                const label faceI = cells[cellI][i];

                maxArea = max
                (
                    maxArea,
                    mesh.isInternalFace(faceI)
                  ? magSf[faceI]
                  : mag(mesh.faceAreas()[faceI])
                );
            }

            cellSize[cellI] = V[cellI]/maxArea;
        }
    }

    labelList patchIDs(patchNames.size());

    forAll(patchNames, i)
    {
        patchIDs[i] = mesh.boundaryMesh().findPatchID(patchNames[i]);

        if (patchIDs[i] < 0)
        {
            FatalErrorIn(args.executable())
                << "Patch " << patchNames[i] << " not found. Valid patches: "
                << mesh.boundaryMesh().names() << exit(FatalError);
        }
    }


    // Null-collision frequency

    const scalar maxGasDensity = gMax(gasDensity.internalField());

    scalar maxSigmaG = 0;

    for (label i = 0; i <= 400; i++)
    {
        const scalar eps = maxCollisionEnergy*Foam::pow(1e-5, 1 - i/400.0);
        const scalar g = Foam::sqrt(2*eps*eCharge/ionMass);

        maxSigmaG = max
        (
            maxSigmaG,
            (sigma.isotropic(eps) + sigma.chargeExchange(eps))*g
        );
    }

    const scalar nuMax = maxGasDensity*maxSigmaG;

    Info<< "Ions produced: " << ionsPerSecond << " 1/s; one tracked ion"
        << " stands for " << weight << " 1/s" << nl
        << "Null-collision frequency " << nuMax << " 1/s" << nl << endl;


    // Scoring

    List<DynamicList<scalar> > scoredEnergy(patchIDs.size());
    List<DynamicList<scalar> > scoredAngle(patchIDs.size());
    List<DynamicList<scalar> > scoredPhase(patchIDs.size());

    scalarField cellTime(mesh.nCells(), 0.0);
    scalarField cellEnergy(mesh.nCells(), 0.0);
    vectorField cellVelocity(mesh.nCells(), vector::zero);

    label nOtherPatch = 0;
    label nTimedOut = 0;
    label nCollisions = 0;
    label nNull = 0;
    label nAboveMax = 0;

    Cloud<ionParticle> cloud(mesh, "ionTrackerCloud", IDLList<ionParticle>());

    const scalar qm = charge/ionMass;
    const scalar dtMax = period/stepsPerPeriod;


    for (label particleI = 0; particleI < nParticles; particleI++)
    {
        // Launch cell, from the cumulative source
        const scalar r = rndGen.scalar01();

        label low = 0;
        label high = mesh.nCells() - 1;

        while (low < high)
        {
            const label mid = (low + high)/2;

            if (cumulative[mid] < r)
            {
                low = mid + 1;
            }
            else
            {
                high = mid;
            }
        }

        const label cell0 = low;

        // Random point in the cell
        vector x0 = mesh.cellCentres()[cell0];
        {
            const boundBox bb
            (
                mesh.cells()[cell0].points(mesh.faces(), mesh.points()),
                false
            );

            for (label attempt = 0; attempt < 100; attempt++)
            {
                vector x = mesh.cellCentres()[cell0];

                for (direction cmpt = 0; cmpt < vector::nComponents; cmpt++)
                {
                    if (directions[cmpt] == 1)
                    {
                        x[cmpt] =
                            bb.min()[cmpt]
                          + (0.001 + 0.998*rndGen.scalar01())
                           *(bb.max()[cmpt] - bb.min()[cmpt]);
                    }
                }

                if (mesh.pointInCell(x, cell0))
                {
                    x0 = x;
                    break;
                }
            }
        }

        // The velocity of a gas atom
        const scalar vThermal =
            Foam::sqrt(kB*gasTemperature[cell0]/gasMass);

        vector U0
        (
            vThermal*rndGen.GaussNormal(),
            vThermal*rndGen.GaussNormal(),
            vThermal*rndGen.GaussNormal()
        );

        ionParticle p(cloud, x0, cell0, U0);
        ionParticle::trackData td(cloud);
        td.keepParticle = true;

        const scalar tLaunch = period*rndGen.scalar01();
        scalar t = tLaunch;

        scalar tFlight = -Foam::log(max(rndGen.scalar01(), VSMALL))/nuMax;

        label stuck = 0;

        while (td.keepParticle)
        {
            if (t - tLaunch > maxTime || stuck > 1000)
            {
                nTimedOut++;
                break;
            }

            const label cellI = p.cell();

            // Field at the position of the ion
            const scalar phase = (t/period - ::floor(t/period))*nPhases;
            const label k0 = min(label(phase), nPhases - 1);
            const label k1 = (k0 + 1) % nPhases;
            const scalar w1 = phase - k0;

            const vector dx0 = p.position() - mesh.cellCentres()[cellI];

            const vector a0 =
                qm
               *(
                    (1 - w1)
                   *(EFields[k0][cellI] + (dx0 & gradEFields[k0][cellI]))
                  + w1
                   *(EFields[k1][cellI] + (dx0 & gradEFields[k1][cellI]))
                );

            vector& U = p.U();

            // Step: to the next collision test, limited by the period, by
            // the size of the cell and by the acceleration
            scalar dt = min(tFlight, dtMax);

            const scalar h = cellFraction*cellSize[cellI];

            dt = min(dt, h/(mag(U) + SMALL));
            dt = min(dt, Foam::sqrt(2*h/(mag(a0) + SMALL)));

            vector Uhalf = U + 0.5*dt*a0;

            vector displacement = Uhalf*dt;

            for (direction cmpt = 0; cmpt < vector::nComponents; cmpt++)
            {
                if (directions[cmpt] != 1)
                {
                    displacement[cmpt] = 0;
                }
            }

            p.stepFraction() = 0;

            const scalar dtDone =
                dt*p.trackToFace(p.position() + displacement, td);

            if (dtDone < 1e-6*dt)
            {
                stuck++;
            }
            else
            {
                stuck = 0;
            }

            t += dtDone;
            tFlight -= dtDone;

            if (!td.keepParticle)
            {
                U += a0*dtDone;
                break;
            }

            // Field at the new position
            const label cellN = p.cell();
            const scalar phaseN = (t/period - ::floor(t/period))*nPhases;
            const label n0 = min(label(phaseN), nPhases - 1);
            const label n1 = (n0 + 1) % nPhases;
            const scalar wn = phaseN - n0;

            const vector dxN = p.position() - mesh.cellCentres()[cellN];

            const vector a1 =
                qm
               *(
                    (1 - wn)
                   *(EFields[n0][cellN] + (dxN & gradEFields[n0][cellN]))
                  + wn
                   *(EFields[n1][cellN] + (dxN & gradEFields[n1][cellN]))
                );

            U += 0.5*(a0 + a1)*dtDone;

            cellTime[cellI] += dtDone;
            cellEnergy[cellI] += 0.5*ionMass*magSqr(U)/eCharge*dtDone;
            cellVelocity[cellI] += U*dtDone;

            if (tFlight <= 0)
            {
                // Collision test with an atom of the gas
                const scalar vt =
                    Foam::sqrt(kB*gasTemperature[cellN]/gasMass);

                const vector Ugas
                (
                    vt*rndGen.GaussNormal(),
                    vt*rndGen.GaussNormal(),
                    vt*rndGen.GaussNormal()
                );

                const scalar g = mag(U - Ugas);
                const scalar eps = 0.5*ionMass*sqr(g)/eCharge;

                const scalar sIso = sigma.isotropic(eps);
                const scalar sCx = sigma.chargeExchange(eps);

                const scalar probability =
                    gasDensity[cellN]*g*(sIso + sCx)/nuMax;

                if (probability > 1)
                {
                    nAboveMax++;
                }

                if (rndGen.scalar01() < probability)
                {
                    nCollisions++;

                    if (rndGen.scalar01()*(sIso + sCx) < sCx)
                    {
                        // Charge exchange: the atom becomes the ion
                        U = Ugas;
                    }
                    else
                    {
                        // Isotropic in the centre-of-mass frame
                        const vector Ucm =
                            (ionMass*U + gasMass*Ugas)/(ionMass + gasMass);

                        U =
                            Ucm
                          + gasMass/(ionMass + gasMass)*g
                           *randomDirection(rndGen);
                    }
                }
                else
                {
                    nNull++;
                }

                tFlight =
                    -Foam::log(max(rndGen.scalar01(), VSMALL))/nuMax;
            }
        }

        if (!td.keepParticle && td.hitPatch >= 0)
        {
            bool scored = false;

            forAll(patchIDs, i)
            {
                if (patchIDs[i] == td.hitPatch)
                {
                    const polyPatch& pp = mesh.boundaryMesh()[td.hitPatch];

                    vector n = pp.faceAreas()[td.hitFace];
                    n /= mag(n);

                    const vector& U = p.U();

                    scoredEnergy[i].append(0.5*ionMass*magSqr(U)/eCharge);

                    scoredAngle[i].append
                    (
                        Foam::acos(min(mag(U & n)/(mag(U) + VSMALL), 1.0))
                       *180/mathematicalConstant::pi
                    );

                    scoredPhase[i].append(t/period - ::floor(t/period));

                    scored = true;
                }
            }

            if (!scored)
            {
                nOtherPatch++;
            }
        }

        if ((particleI + 1) % max(nParticles/10, 1) == 0)
        {
            Info<< "    " << particleI + 1 << " ions tracked" << endl;
        }
    }


    // Output

    const fileName outputDir(runTime.path()/"ionTracker");
    mkDir(outputDir);

    OFstream summary(outputDir/"summary.dat");

    summary
        << "# ionTracker: ion " << ionName << ", " << nParticles
        << " ions tracked, fields of the period from " << startTime << " s"
        << nl
        << "# collisions " << nCollisions << ", null collisions " << nNull
        << ", above the null-collision frequency " << nAboveMax << nl
        << "# reached other patches " << nOtherPatch << ", not finished "
        << nTimedOut << nl
        << "# patch  ions  fraction  flux[1/m2/s]  meanEnergy[eV]"
        << "  meanAngle[deg]" << nl;

    Info<< nl << "Collisions " << nCollisions << ", null collisions " << nNull
        << ", tests above the null-collision frequency " << nAboveMax << nl
        << "Ions that reached other patches " << nOtherPatch
        << ", not finished after maxTime " << nTimedOut << nl << endl;

    forAll(patchIDs, i)
    {
        const scalarList& energy = scoredEnergy[i];
        const scalarList& angle = scoredAngle[i];
        const scalarList& phase = scoredPhase[i];

        const label n = energy.size();

        scalar meanEnergy = 0;
        scalar meanAngle = 0;
        scalar largest = 0;

        forAll(energy, j)
        {
            meanEnergy += energy[j];
            meanAngle += angle[j];
            largest = max(largest, energy[j]);
        }

        if (n > 0)
        {
            meanEnergy /= n;
            meanAngle /= n;
        }

        const scalar area = gSum(mesh.magSf().boundaryField()[patchIDs[i]]);
        const scalar flux = weight*n/area;

        Info<< patchNames[i] << ": " << n << " ions, flux " << flux
            << " 1/m2/s, mean energy " << meanEnergy << " eV, mean angle "
            << meanAngle << " deg, largest energy " << largest << " eV"
            << endl;

        summary
            << patchNames[i] << tab << n << tab << scalar(n)/nParticles
            << tab << flux << tab << meanEnergy << tab << meanAngle << nl;

        if (n == 0)
        {
            continue;
        }

        const scalar eMax = (maxEnergy > 0 ? maxEnergy : 1.0001*largest);
        const scalar dE = eMax/nEnergyBins;
        const scalar dA = 90.0/nAngleBins;

        scalarField fE(nEnergyBins, 0.0);
        scalarField fA(nAngleBins, 0.0);
        List<scalarField> fEA(nEnergyBins, scalarField(nAngleBins, 0.0));

        forAll(energy, j)
        {
            const label iE = min(label(energy[j]/dE), nEnergyBins - 1);
            const label iA = min(label(angle[j]/dA), nAngleBins - 1);

            fE[iE] += 1;
            fA[iA] += 1;
            fEA[iE][iA] += 1;
        }

        {
            OFstream os(outputDir/(patchNames[i] + "_energy.dat"));

            os  << "# energy [eV]" << tab << "f(E) [1/eV]" << nl;

            forAll(fE, iE)
            {
                os  << (iE + 0.5)*dE << tab << fE[iE]/(n*dE) << nl;
            }
        }

        {
            OFstream os(outputDir/(patchNames[i] + "_angle.dat"));

            os  << "# angle to the normal [deg]" << tab << "f(angle) [1/deg]"
                << nl;

            forAll(fA, iA)
            {
                os  << (iA + 0.5)*dA << tab << fA[iA]/(n*dA) << nl;
            }
        }

        {
            OFstream os(outputDir/(patchNames[i] + "_energyAngle.dat"));

            os  << "# rows: energy bins of " << dE << " eV from 0;"
                << " columns: angle bins of " << dA << " deg from 0;"
                << " f(E, angle) [1/eV/deg]" << nl;

            forAll(fEA, iE)
            {
                forAll(fEA[iE], iA)
                {
                    os  << fEA[iE][iA]/(n*dE*dA) << ' ';
                }

                os  << nl;
            }
        }

        if (writeParticles)
        {
            OFstream os(outputDir/(patchNames[i] + "_particles.dat"));

            os  << "# energy [eV]" << tab << "angle [deg]" << tab
                << "phase of arrival" << nl;

            forAll(energy, j)
            {
                os  << energy[j] << tab << angle[j] << tab << phase[j] << nl;
            }
        }
    }


    // Fields of the tracked ions

    volScalarField ionTrackerDensity
    (
        IOobject
        (
            "ionTrackerDensity",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    );

    volScalarField ionTrackerEnergy
    (
        IOobject
        (
            "ionTrackerEnergy",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    );

    volVectorField ionTrackerVelocity
    (
        IOobject
        (
            "ionTrackerVelocity",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimless, vector::zero),
        zeroGradientFvPatchVectorField::typeName
    );

    forAll(cellTime, cellI)
    {
        ionTrackerDensity[cellI] = weight*cellTime[cellI]/V[cellI];

        if (cellTime[cellI] > 0)
        {
            ionTrackerEnergy[cellI] = cellEnergy[cellI]/cellTime[cellI];
            ionTrackerVelocity[cellI] = cellVelocity[cellI]/cellTime[cellI];
        }
    }

    ionTrackerDensity.correctBoundaryConditions();
    ionTrackerEnergy.correctBoundaryConditions();
    ionTrackerVelocity.correctBoundaryConditions();

    ionTrackerDensity.write();
    ionTrackerEnergy.write();
    ionTrackerVelocity.write();

    Info<< nl << "Distributions written to " << outputDir << nl
        << "Fields ionTrackerDensity, ionTrackerEnergy and ionTrackerVelocity"
        << " written to time " << runTime.timeName() << nl << nl
        << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
