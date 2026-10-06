/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

Application
    electronBoltzmann

Description
    Electron transport and rate coefficients from cross sections: the
    Boltzmann equation for the electrons in a uniform field, in the two-term
    approximation (Hagelaar and Pitchford, Plasma Sources Sci. Technol. 14
    (2005) 722), for a list of reduced fields E/N.

    With eps the electron energy [eV], F0 the isotropic part of the energy
    distribution, normalised to int sqrt(eps) F0 d(eps) = 1, and
    gamma = sqrt(2e/m), the equation solved is

        d/d(eps) [ W F0 - D dF0/d(eps) ] = S

        W = -gamma eps^2 sum_k x_k (2m/M_k) sigma_k,elastic
        D = gamma/3 (E/N)^2 eps/sigma_m
          + gamma kB T/e eps^2 sum_k x_k (2m/M_k) sigma_k,elastic
        S = gamma sum_inelastic x_k [ (eps + u) sigma(eps + u) F0(eps + u)
                                    - eps sigma(eps) F0(eps) ]

    with sigma_m the total momentum-transfer cross section of the mixture
    and u the energy loss of each inelastic process. It is discretised on a
    uniform energy grid with the exponential (Scharfetter-Gummel) flux and
    solved directly; the upper end of the grid is adapted to the
    distribution. Ionisation is treated as an energy loss (the electron
    that is released is not added) and attachment only gives a rate
    coefficient, so the result is for conditions where the electron number
    changes slowly on the scale of the collisions; there are no
    electron-electron collisions.

    With "method monteCarlo;" the coefficients are instead from a Monte
    Carlo simulation of a swarm of electrons in the same uniform field
    (null-collision method), which does not use the two-term approximation
    and adds the electrons released by ionisation and removes the attached
    ones. Scattering is isotropic, with the elastic (momentum-transfer)
    cross section of the file as the elastic cross section. The swarm
    starts from the two-term distribution. The drift velocity is the mean
    velocity of the electrons and the diffusion coefficient the transverse
    one. transportErrors.dat holds the standard errors and
    transportTwoTerm.dat the two-term results.

    Cross sections are read from a file in the LXCat format: blocks ELASTIC
    or EFFECTIVE (with the mass ratio m/M), EXCITATION and IONIZATION (with
    the energy loss) and ATTACHMENT, for the targets listed.

    system/boltzmannDict:

    \verbatim
    file            "constant/Ar.txt";      // LXCat format
    gases           ((Ar 1.0));             // (target fraction) ...
    gasTemperature  300;                    // [K]

    reducedField    { min 0.01; max 1000; n 61; }   // [Td], logarithmic
    // reducedField (1 10 100);                     // or a list

    nEnergyCells    400;
    outputDir       "constant/boltzmann";

    method          twoTerm;                // or monteCarlo, with
    // nElectrons      1000;
    // sampleTimes     20;      // sampling time in energy relaxation times
    // equilibrationTimes 2;    // time before sampling
    // energySharing   0.5;     // share of the electron released in an
    //                          // ionisation in the energy that is left
    // seed            1234;

    // Tables of the fluid solver to complete (optional): the points that
    // are there are kept, and computed points are added below the first
    // and above the last electron temperature of each table
    // extend
    // {
    //     mobility    "constant/mu_electron";
    //     diffusion   "constant/D_electron";
    //     rates       (("constant/reaction_2" 3));   // (file process)
    // }
    \endverbatim

    Output in outputDir: transport.dat with, for each E/N, the mean energy,
    the electron temperature Te = 2/3 mean energy, mu N, D N, the energy
    mobility and diffusion coefficient times N, and the rate coefficients of
    all processes [m3/s]; and tables against Te in the formats of the fluid
    solver: mu_electron and D_electron (Te [K], mu N or D N) and rate_<n>
    (Te [K], ln of the rate coefficient in m3/kmol/s) for process n in the
    order of the file. distribution_<E/N>.dat holds the distributions.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "foamTime.H"
#include "IOdictionary.H"
#include "IFstream.H"
#include "OFstream.H"
#include "IStringStream.H"
#include "OStringStream.H"
#include "Tuple2.H"
#include "DynamicList.H"
#include "scalarField.H"
#include "scalarSquareMatrix.H"
#include "OSspecific.H"
#include "Random.H"
#include "mathematicalConstants.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- A collision process of the cross section file
struct process
{
    word type;          // ELASTIC, EFFECTIVE, EXCITATION, IONIZATION, ATTACHMENT
    string target;      // line after the keyword
    scalar parameter;   // m/M, or the energy loss [eV]
    scalar fraction;    // mole fraction of the target
    scalarField energy; // [eV]
    scalarField sigma;  // [m2]

    //- Linear interpolation; zero below the first point of a process with
    //  a threshold, constant beyond the ends otherwise
    scalar operator()(const scalar eps) const
    {
        const label n = energy.size();

        if (eps <= energy[0])
        {
            return (type == "ELASTIC" || type == "EFFECTIVE") ? sigma[0] : 0;
        }

        if (eps >= energy[n - 1])
        {
            return sigma[n - 1];
        }

        label low = 0;
        label high = n - 1;

        while (high - low > 1)
        {
            const label mid = (low + high)/2;

            if (energy[mid] <= eps)
            {
                low = mid;
            }
            else
            {
                high = mid;
            }
        }

        const scalar w = (eps - energy[low])/(energy[high] - energy[low]);

        return (1 - w)*sigma[low] + w*sigma[high];
    }
};


string trim(const string& line)
{
    const string::size_type first = line.find_first_not_of(" \t\r");

    if (first == string::npos)
    {
        return string("");
    }

    const string::size_type last = line.find_last_not_of(" \t\r");

    return string(line.substr(first, last - first + 1));
}


bool isDashes(const string& s)
{
    return s.size() >= 5 && s.substr(0, 5) == "-----";
}


//- Read the processes of the listed targets from a file in the LXCat format
void readLXCat
(
    const fileName& file,
    const List<Tuple2<word, scalar> >& gases,
    PtrList<process>& processes
)
{
    IFstream is(file);

    if (!is.good())
    {
        FatalErrorIn("readLXCat") << "Cannot read " << file
            << exit(FatalError);
    }

    DynamicList<process*> found;

    while (is.good())
    {
        string line;
        is.getLine(line);

        const string keyword(trim(line));

        if
        (
            keyword != "ELASTIC" && keyword != "EFFECTIVE"
         && keyword != "EXCITATION" && keyword != "IONIZATION"
         && keyword != "ATTACHMENT" && keyword != "ROTATION"
        )
        {
            continue;
        }

        string targetLine;
        is.getLine(targetLine);
        targetLine = trim(targetLine);

        // The target is the first word of the line
        const string target
        (
            targetLine.substr(0, targetLine.find_first_of(" \t"))
        );

        scalar fraction = -1;

        forAll(gases, i)
        {
            if (target == gases[i].first())
            {
                fraction = gases[i].second();
            }
        }

        scalar parameter = 0;

        if (keyword != "ATTACHMENT")
        {
            string parameterLine;
            is.getLine(parameterLine);

            IStringStream ps(trim(parameterLine));
            ps >> parameter;
        }

        // The table between the two lines of dashes
        DynamicList<scalar> e;
        DynamicList<scalar> s;

        bool inTable = false;

        while (is.good())
        {
            string row;
            is.getLine(row);
            row = trim(row);

            if (isDashes(row))
            {
                if (inTable)
                {
                    break;
                }

                inTable = true;
            }
            else if (inTable && row.size())
            {
                IStringStream rs(row);
                scalar a, b;
                rs >> a >> b;

                // Energies must increase
                if (e.empty() || a > e[e.size() - 1])
                {
                    e.append(a);
                    s.append(b);
                }
            }
        }

        if (fraction <= 0 || e.size() < 2)
        {
            continue;
        }

        process* pPtr = new process;

        pPtr->type = (keyword == "ROTATION" ? word("EXCITATION") : word(keyword));
        pPtr->target = targetLine;
        pPtr->parameter = parameter;
        pPtr->fraction = fraction;
        pPtr->energy = scalarField(e);
        pPtr->sigma = scalarField(s);

        found.append(pPtr);
    }

    processes.setSize(found.size());

    forAll(found, i)
    {
        processes.set(i, found[i]);
    }
}


//- Bernoulli function x/(exp(x) - 1)
scalar bernoulli(const scalar x)
{
    if (mag(x) < 1e-6)
    {
        return 1 - 0.5*x;
    }
    else if (x > 500)
    {
        return 0;
    }
    else if (x < -500)
    {
        return -x;
    }

    return x/(Foam::exp(x) - 1);
}


//- Momentum-transfer cross section of the mixture and the elastic energy
//  exchange coefficient sum_k x_k (2m/M_k) sigma_k,elastic at an energy
void mixtureCrossSections
(
    const PtrList<process>& processes,
    const List<Tuple2<word, scalar> >& gases,
    const scalar eps,
    scalar& sigmaM,
    scalar& sigmaEps
)
{
    sigmaM = 0;
    sigmaEps = 0;

    forAll(gases, g)
    {
        scalar elastic = 0;
        scalar effective = -1;
        scalar inelastic = 0;
        scalar massRatio = 0;

        forAll(processes, i)
        {
            const process& p = processes[i];

            const string target
            (
                p.target.substr(0, p.target.find_first_of(" \t"))
            );

            if (target != gases[g].first())
            {
                continue;
            }

            if (p.type == "ELASTIC")
            {
                elastic = p(eps);
                massRatio = p.parameter;
            }
            else if (p.type == "EFFECTIVE")
            {
                effective = p(eps);
                massRatio = p.parameter;
            }
            else if (p.type == "EXCITATION" || p.type == "IONIZATION")
            {
                inelastic += p(eps);
            }
        }

        const scalar x = gases[g].second();

        if (effective >= 0)
        {
            // The effective cross section includes the inelastic ones
            sigmaM += x*effective;
            sigmaEps += x*2*massRatio*max(effective - inelastic, 0.0);
        }
        else
        {
            sigmaM += x*(elastic + inelastic);
            sigmaEps += x*2*massRatio*elastic;
        }
    }
}


//- Solve for the distribution on a uniform grid of n cells up to epsMax
void solveDistribution
(
    const PtrList<process>& processes,
    const List<Tuple2<word, scalar> >& gases,
    const scalar EN,
    const scalar kTe,
    const scalar epsMax,
    scalarField& F
)
{
    const scalar gamma = Foam::sqrt(2*1.602176634e-19/9.1093837015e-31);

    const label n = F.size();
    const scalar dEps = epsMax/n;

    scalarSquareMatrix A(n, 0.0);

    // Fluxes through the faces between the cells, W F - D dF/d(eps)
    for (label f = 1; f < n; f++)
    {
        const scalar eps = f*dEps;

        scalar sigmaM, sigmaEps;
        mixtureCrossSections(processes, gases, eps, sigmaM, sigmaEps);

        const scalar W = -gamma*sqr(eps)*sigmaEps;

        const scalar D =
            gamma/3*sqr(EN)*eps/max(sigmaM, 1e-30)
          + gamma*kTe*sqr(eps)*sigmaEps;

        const scalar z = W*dEps/max(D, 1e-300);

        // Flux from cell f - 1 to cell f
        const scalar aLow = D/dEps*bernoulli(-z);
        const scalar aHigh = D/dEps*bernoulli(z);

        // d(flux)/d(eps) = S: flux out of f - 1, into f
        A[f - 1][f - 1] += aLow;
        A[f - 1][f] -= aHigh;
        A[f][f - 1] -= aLow;
        A[f][f] += aHigh;
    }

    // Inelastic collisions: out of cell j, into the cell at eps_j - u
    forAll(processes, i)
    {
        const process& p = processes[i];

        if (p.type != "EXCITATION" && p.type != "IONIZATION")
        {
            continue;
        }

        const scalar u = p.parameter;

        for (label j = 0; j < n; j++)
        {
            const scalar eps = (j + 0.5)*dEps;

            if (eps <= u)
            {
                continue;
            }

            const scalar rate = gamma*p.fraction*eps*p(eps)*dEps;

            const label target = min(max(label((eps - u)/dEps), 0), n - 1);

            if (target != j)
            {
                A[j][j] += rate;
                A[target][j] -= rate;
            }
        }
    }

    // The equations are linearly dependent (electrons are conserved): the
    // first is replaced by the normalisation
    scalarField b(n, 0.0);

    for (label j = 0; j < n; j++)
    {
        A[0][j] = Foam::sqrt((j + 0.5)*dEps)*dEps;
    }

    b[0] = 1;

    // Scale the rows for the elimination
    for (label i = 1; i < n; i++)
    {
        scalar rowMax = 0;

        for (label j = 0; j < n; j++)
        {
            rowMax = max(rowMax, mag(A[i][j]));
        }

        if (rowMax > 0)
        {
            for (label j = 0; j < n; j++)
            {
                A[i][j] /= rowMax;
            }
        }
        else
        {
            A[i][i] = 1;
        }
    }

    scalarSquareMatrix::LUsolve(A, b);

    F = b;
}


//- Results for one reduced field
struct swarm
{
    scalar EN;              // [V m2]
    scalar meanEnergy;      // [eV]
    scalar muN;             // [1/(V m s)]
    scalar DN;              // [1/(m s)]
    scalar muEpsN;
    scalar DEpsN;
    scalarField rates;      // [m3/s]
    scalar epsMax;
};


//- Write the table of the results against E/N
void writeTransport
(
    const fileName& file,
    const List<swarm>& results,
    const label nProcesses,
    const string& note
)
{
    const scalar Td = 1e-21;

    OFstream os(file);

    os  << "# " << note.c_str() << nl
        << "# E/N [Td]" << tab << "mean energy [eV]" << tab << "Te [K]"
        << tab << "mu N [1/(V m s)]" << tab << "D N [1/(m s)]" << tab
        << "muEps N" << tab << "DEps N";

    for (label i = 0; i < nProcesses; i++)
    {
        os  << tab << "k" << i << " [m3/s]";
    }

    os  << nl;

    forAll(results, k)
    {
        const swarm& r = results[k];

        os  << r.EN/Td << tab << r.meanEnergy << tab
            << 2.0/3.0*r.meanEnergy*1.602176634e-19/1.380649e-23 << tab
            << r.muN << tab << r.DN << tab << r.muEpsN << tab << r.DEpsN;

        forAll(r.rates, i)
        {
            os  << tab << r.rates[i];
        }

        os  << nl;
    }
}


//- Settings of the Monte Carlo method
struct monteCarloSettings
{
    label nElectrons;
    scalar gasDensity;
    scalar energySharing;
    scalar equilibrationTimes;
    scalar sampleTimes;
    scalar minCollisions;
    label nBlocks;
    label nSpeedCells;
};


//- An electron of the swarm
struct swarmElectron
{
    vector v;
    scalar x;
    scalar y;
    scalar t;
    bool alive;
};


vector isotropicDirection(Random& rndGen)
{
    const scalar cosTheta = 2*rndGen.scalar01() - 1;
    const scalar sinTheta = Foam::sqrt(max(1 - sqr(cosTheta), 0.0));
    const scalar phi = 2*mathematicalConstant::pi*rndGen.scalar01();

    return vector(sinTheta*Foam::cos(phi), sinTheta*Foam::sin(phi), cosTheta);
}


//- Monte Carlo simulation of a swarm of electrons in a uniform field and
//  gas (null-collision method). The electrons start from the two-term
//  distribution F (cells of width dEpsF); r holds the two-term results on
//  entry and the Monte Carlo results on return, err their standard errors
//  from the scatter of the blocks of the sampling time.
//
//  Scattering is isotropic. The elastic cross section of the file (a
//  momentum-transfer cross section) is used as the elastic cross section,
//  and the gas atoms have a Maxwellian velocity in elastic collisions. An
//  excitation takes its energy loss from the electron. An ionisation
//  releases an electron: the energy that is left is shared between the two
//  (energySharing is the fraction of the released one). An attachment
//  removes the electron. The number of simulated electrons is kept between
//  half and twice nElectrons.
//
//  The drift velocity is the mean velocity of the electrons (flux drift
//  velocity) and the diffusion coefficient is the transverse one, from the
//  spread of the swarm across the field; the rate coefficients are
//  averages of sigma v over the electrons.
void monteCarlo
(
    const PtrList<process>& processes,
    const List<Tuple2<word, scalar> >& gases,
    const scalar gasTemperature,
    const scalarField& F,
    const scalar dEpsF,
    const monteCarloSettings& mc,
    Random& rndGen,
    swarm& r,
    swarm& err,
    scalarField& distribution,
    scalar& speedCell
)
{
    const scalar eCharge = 1.602176634e-19;
    const scalar eMass = 9.1093837015e-31;
    const scalar kB = 1.380649e-23;
    const scalar gamma = Foam::sqrt(2*eCharge/eMass);

    const scalar EN = r.EN;
    const scalar N = mc.gasDensity;
    const scalar aZ = eCharge*EN*N/eMass;

    // Speed grid; the cross sections are constant in its cells
    const label nGrid = mc.nSpeedCells;
    const scalar vMax = gamma*Foam::sqrt(max(1.5*F.size()*dEpsF, 1.0));
    const scalar dV = vMax/nGrid;

    speedCell = dV;

    // Collision processes of the simulation: 0 elastic, 1 excitation,
    // 2 ionisation, 3 attachment
    DynamicList<label> kind;
    DynamicList<scalar> loss;
    DynamicList<scalar> massRatio;
    DynamicList<scalarField*> sigmaPtrs;

    forAll(gases, g)
    {
        // Elastic: the effective cross section less the inelastic ones
        scalarField elastic(nGrid, 0.0);
        scalarField inelastic(nGrid, 0.0);
        bool effective = false;
        scalar ratio = 0;

        forAll(processes, i)
        {
            const process& p = processes[i];

            const string target
            (
                p.target.substr(0, p.target.find_first_of(" \t"))
            );

            if (target != gases[g].first())
            {
                continue;
            }

            scalarField* sPtr = new scalarField(nGrid, 0.0);

            for (label j = 0; j < nGrid; j++)
            {
                (*sPtr)[j] = p.fraction*p(sqr((j + 0.5)*dV/gamma));
            }

            if (p.type == "ELASTIC" || p.type == "EFFECTIVE")
            {
                elastic = *sPtr;
                effective = (p.type == "EFFECTIVE");
                ratio = p.parameter;
                delete sPtr;
            }
            else
            {
                if (p.type != "ATTACHMENT")
                {
                    inelastic += *sPtr;
                }

                kind.append
                (
                    p.type == "EXCITATION" ? 1
                  : p.type == "IONIZATION" ? 2
                  : 3
                );
                loss.append(p.parameter);
                massRatio.append(0);
                sigmaPtrs.append(sPtr);
            }
        }

        if (effective)
        {
            elastic = max(elastic - inelastic, scalar(0));
        }

        kind.append(0);
        loss.append(0);
        massRatio.append(ratio);
        sigmaPtrs.append(new scalarField(elastic));
    }

    const label nProc = kind.size();

    PtrList<scalarField> sigma(nProc);

    forAll(sigma, k)
    {
        sigma.set(k, sigmaPtrs[k]);
    }

    // Largest collision frequency up to each cell
    scalarField sigmaTotal(nGrid, 0.0);

    forAll(sigma, k)
    {
        sigmaTotal += sigma[k];
    }

    scalarField nuMaxUpTo(nGrid, 0.0);

    forAll(nuMaxUpTo, j)
    {
        nuMaxUpTo[j] = max(N*sigmaTotal[j]*(j + 1)*dV, SMALL);

        if (j > 0)
        {
            nuMaxUpTo[j] = max(nuMaxUpTo[j], nuMaxUpTo[j - 1]);
        }
    }

    const scalar vFloor = gamma*Foam::sqrt(0.05);

    // Time scales from the two-term results: energy relaxation by the
    // field and by elastic collisions
    scalar sigmaM, sigmaEps;
    mixtureCrossSections(processes, gases, r.meanEnergy, sigmaM, sigmaEps);

    scalar tau = r.meanEnergy/max(sqr(EN)*r.muN*N, SMALL);

    if (sigmaEps > 0)
    {
        tau = min(tau, 1/(N*gamma*Foam::sqrt(r.meanEnergy)*sigmaEps));
    }

    const scalar nuMean =
        N*gamma*Foam::sqrt(r.meanEnergy)*max(sigmaM, 1e-30);

    const scalar tSample =
        max(mc.sampleTimes*tau, mc.minCollisions/nuMean);

    const scalar tBlock = tSample/mc.nBlocks;

    const label nEquilibration =
        max(label(mc.equilibrationTimes*tau/tBlock + 0.999), 1);

    // Electrons from the two-term distribution
    DynamicList<swarmElectron> electrons(2*mc.nElectrons);

    {
        scalarField cumulative(F.size());
        scalar sum = 0;

        forAll(F, j)
        {
            sum += Foam::sqrt((j + 0.5)*dEpsF)*max(F[j], 0.0);
            cumulative[j] = sum;
        }

        for (label n = 0; n < mc.nElectrons; n++)
        {
            const scalar u = rndGen.scalar01()*sum;

            label low = 0;
            label high = F.size() - 1;

            while (low < high)
            {
                const label mid = (low + high)/2;

                if (cumulative[mid] < u)
                {
                    low = mid + 1;
                }
                else
                {
                    high = mid;
                }
            }

            swarmElectron e;

            e.v =
                gamma*Foam::sqrt((low + rndGen.scalar01())*dEpsF)
               *isotropicDirection(rndGen);
            e.x = 0;
            e.y = 0;
            e.t = 0;
            e.alive = true;

            electrons.append(e);
        }
    }

    // Sums of the blocks
    scalarField blockTime(mc.nBlocks, 0.0);
    scalarField blockEnergy(mc.nBlocks, 0.0);
    scalarField blockDz(mc.nBlocks, 0.0);
    scalarField blockD(mc.nBlocks, 0.0);

    // Time spent in each speed cell, and its integral with the speed
    List<scalarField> blockH(mc.nBlocks, scalarField(nGrid, 0.0));
    List<scalarField> blockHv(mc.nBlocks, scalarField(nGrid, 0.0));

    scalar nCollisions = 0;
    scalar nTests = 0;
    scalar nBeyond = 0;

    for (label block = -nEquilibration; block < mc.nBlocks; block++)
    {
        const bool sampling = (block >= 0);

        forAll(electrons, n)
        {
            electrons[n].t = 0;
            electrons[n].x = 0;
            electrons[n].y = 0;
        }

        scalar time = 0;
        scalar energyTime = 0;
        scalar dz = 0;

        // The list grows when electrons are released
        for (label n = 0; n < electrons.size(); n++)
        {
            vector v = electrons[n].v;
            scalar x = electrons[n].x;
            scalar y = electrons[n].y;
            scalar t = electrons[n].t;
            bool alive = true;

            while (alive)
            {
                const scalar speed = mag(v);

                // The collision frequency is bounded while the speed
                // stays below the upper edge of a cell above it
                const label cap =
                    min(label(max(2*speed, vFloor)/dV), nGrid - 1);

                const scalar nuMax = nuMaxUpTo[cap];

                const scalar tLimit =
                (
                    cap == nGrid - 1
                  ? GREAT
                  : ((cap + 1)*dV - speed)/aZ
                );

                const scalar tFree =
                    -Foam::log(max(rndGen.scalar01(), 1e-300))/nuMax;

                const scalar tLeft = tBlock - t;

                const scalar dt = min(tFree, min(tLimit, tLeft));

                // Free flight
                if (sampling)
                {
                    time += dt;

                    energyTime +=
                        0.5*eMass/eCharge
                       *(
                            magSqr(v)*dt + v.z()*aZ*sqr(dt)
                          + sqr(aZ)*pow3(dt)/3
                        );

                    dz += v.z()*dt + 0.5*aZ*sqr(dt);
                }

                x += v.x()*dt;
                y += v.y()*dt;
                v.z() += aZ*dt;
                t += dt;

                if (tLeft <= tFree && tLeft <= tLimit)
                {
                    break;
                }
                else if (tLimit < tFree)
                {
                    continue;
                }

                // Test for a collision
                const scalar speedNow = mag(v);
                label j = label(speedNow/dV);

                if (j >= nGrid)
                {
                    j = nGrid - 1;
                    nBeyond++;
                }

                if (sampling)
                {
                    blockH[block][j] += 1/nuMax;
                    blockHv[block][j] += speedNow/nuMax;
                }

                nTests++;

                const scalar u = rndGen.scalar01()*nuMax/(N*speedNow + VSMALL);

                scalar sum = 0;
                label chosen = -1;

                for (label k = 0; k < nProc; k++)
                {
                    sum += sigma[k][j];

                    if (u < sum)
                    {
                        chosen = k;
                        break;
                    }
                }

                if (chosen < 0)
                {
                    continue;
                }

                const scalar eps = sqr(speedNow/gamma);

                if (kind[chosen] == 0)
                {
                    // Elastic collision with an atom of the gas
                    const scalar ratio = massRatio[chosen];

                    if (ratio > 0)
                    {
                        const scalar atomMass = eMass/ratio;

                        vector V(vector::zero);

                        if (gasTemperature > 0)
                        {
                            const scalar vt =
                                Foam::sqrt(kB*gasTemperature/atomMass);

                            V = vt*vector
                            (
                                rndGen.GaussNormal(),
                                rndGen.GaussNormal(),
                                rndGen.GaussNormal()
                            );
                        }

                        const vector vCentre = (ratio*v + V)/(1 + ratio);

                        v =
                            vCentre
                          + mag(v - V)/(1 + ratio)
                           *isotropicDirection(rndGen);
                    }
                    else
                    {
                        v = speedNow*isotropicDirection(rndGen);
                    }
                }
                else if (kind[chosen] == 1)
                {
                    if (eps <= loss[chosen])
                    {
                        continue;
                    }

                    v =
                        gamma*Foam::sqrt(eps - loss[chosen])
                       *isotropicDirection(rndGen);
                }
                else if (kind[chosen] == 2)
                {
                    if (eps <= loss[chosen])
                    {
                        continue;
                    }

                    const scalar left = eps - loss[chosen];

                    swarmElectron released;

                    released.v =
                        gamma*Foam::sqrt(mc.energySharing*left)
                       *isotropicDirection(rndGen);
                    released.x = x;
                    released.y = y;
                    released.t = t;
                    released.alive = true;

                    electrons.append(released);

                    v =
                        gamma*Foam::sqrt((1 - mc.energySharing)*left)
                       *isotropicDirection(rndGen);
                }
                else
                {
                    alive = false;
                }

                nCollisions++;
            }

            electrons[n].v = v;
            electrons[n].x = x;
            electrons[n].y = y;
            electrons[n].t = t;
            electrons[n].alive = alive;
        }

        // Remove the attached electrons
        {
            label kept = 0;

            forAll(electrons, n)
            {
                if (electrons[n].alive)
                {
                    electrons[kept++] = electrons[n];
                }
            }

            electrons.setSize(kept);
        }

        if (electrons.empty())
        {
            FatalErrorIn("monteCarlo")
                << "All the electrons are attached at E/N " << EN/1e-21
                << " Td" << exit(FatalError);
        }

        if (sampling)
        {
            scalar spread = 0;

            forAll(electrons, n)
            {
                spread += sqr(electrons[n].x) + sqr(electrons[n].y);
            }

            blockTime[block] = time;
            blockEnergy[block] = energyTime;
            blockDz[block] = dz;
            blockD[block] = spread/electrons.size()/(4*tBlock);
        }

        // Keep the number of electrons between half and twice the nominal
        while (electrons.size() > 2*mc.nElectrons)
        {
            const label nOld = electrons.size();

            for (label n = 0; n < nOld/2; n++)
            {
                const label other = rndGen.integer(n, nOld - 1);

                const swarmElectron e = electrons[n];
                electrons[n] = electrons[other];
                electrons[other] = e;
            }

            electrons.setSize(nOld/2);
        }

        while (electrons.size() < mc.nElectrons/2)
        {
            const label nOld = electrons.size();

            for (label n = 0; n < nOld; n++)
            {
                electrons.append(electrons[n]);
            }
        }
    }

    // Means and standard errors over the blocks
    const label nB = mc.nBlocks;
    const label nOut = processes.size();

    // Cross sections of the processes of the file in the speed cells
    List<scalarField> sigmaOut(nOut, scalarField(nGrid, 0.0));

    forAll(processes, i)
    {
        for (label j = 0; j < nGrid; j++)
        {
            sigmaOut[i][j] = processes[i](sqr((j + 0.5)*dV/gamma));
        }
    }

    scalarField energy(nB);
    scalarField mobility(nB);
    scalarField diffusion(nB);
    List<scalarField> rates(nOut, scalarField(nB, 0.0));

    for (label b = 0; b < nB; b++)
    {
        energy[b] = blockEnergy[b]/blockTime[b];
        mobility[b] = blockDz[b]/blockTime[b]/EN;
        diffusion[b] = blockD[b]*N;

        const scalar weight = sum(blockH[b]);

        forAll(processes, i)
        {
            rates[i][b] = sum(sigmaOut[i]*blockHv[b])/weight;
        }
    }

    r.meanEnergy = average(energy);
    r.muN = average(mobility);
    r.DN = average(diffusion);

    err = r;

    const scalar root = Foam::sqrt(scalar(max(nB*(nB - 1), 1)));

    err.meanEnergy = Foam::sqrt(sum(sqr(energy - r.meanEnergy)))/root;
    err.muN = Foam::sqrt(sum(sqr(mobility - r.muN)))/root;
    err.DN = Foam::sqrt(sum(sqr(diffusion - r.DN)))/root;

    forAll(processes, i)
    {
        r.rates[i] = average(rates[i]);
        err.rates[i] = Foam::sqrt(sum(sqr(rates[i] - r.rates[i])))/root;
    }

    // Distribution F0 in the speed cells
    distribution.setSize(nGrid);
    distribution = 0.0;

    for (label b = 0; b < nB; b++)
    {
        distribution += blockH[b];
    }

    distribution /= sum(distribution);

    forAll(distribution, j)
    {
        const scalar eps = sqr((j + 0.5)*dV/gamma);
        const scalar dEps = sqr((j + 1)*dV/gamma) - sqr(j*dV/gamma);

        distribution[j] /= Foam::sqrt(eps)*dEps;
    }

    Info<< "    Monte Carlo: mean energy " << r.meanEnergy << " +- "
        << err.meanEnergy << " eV, mu N " << r.muN << " +- " << err.muN
        << ", D N " << r.DN << " +- " << err.DN << "; " << nCollisions
        << " collisions, " << nTests - nCollisions << " null, "
        << electrons.size() << " electrons at the end";

    if (nBeyond > 0)
    {
        Info<< "; " << nBeyond << " tests beyond the speed grid";
    }

    Info<< endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"

    const scalar eCharge = 1.602176634e-19;
    const scalar kB = 1.380649e-23;
    const scalar gamma = Foam::sqrt(2*eCharge/9.1093837015e-31);
    const scalar Td = 1e-21;

    IOdictionary dict
    (
        IOobject
        (
            "boltzmannDict",
            runTime.system(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    fileName file(dict.lookup("file"));
    file.expand();

    if (file[0] != '/')
    {
        file = runTime.path()/file;
    }

    const List<Tuple2<word, scalar> > gases(dict.lookup("gases"));
    const scalar gasTemperature = readScalar(dict.lookup("gasTemperature"));
    const label nCells = dict.lookupOrDefault<label>("nEnergyCells", 400);

    fileName outputDir
    (
        dict.lookupOrDefault<fileName>("outputDir", "constant/boltzmann")
    );

    if (outputDir[0] != '/')
    {
        outputDir = runTime.path()/outputDir;
    }

    // Reduced fields [Td]
    scalarList fields;

    if (dict.isDict("reducedField"))
    {
        const dictionary& rf = dict.subDict("reducedField");

        const scalar lo = readScalar(rf.lookup("min"));
        const scalar hi = readScalar(rf.lookup("max"));
        const label n = readLabel(rf.lookup("n"));

        fields.setSize(n);

        forAll(fields, i)
        {
            fields[i] =
                (n == 1 ? lo : lo*Foam::pow(hi/lo, scalar(i)/(n - 1)));
        }
    }
    else
    {
        fields = scalarList(dict.lookup("reducedField"));
    }

    PtrList<process> processes;
    readLXCat(file, gases, processes);

    Info<< "Cross sections from " << file << nl;

    forAll(processes, i)
    {
        Info<< "    " << i << "  " << processes[i].type << "  "
            << processes[i].target.c_str() << "  parameter "
            << processes[i].parameter << ", " << processes[i].energy.size()
            << " points" << nl;
    }

    if (processes.empty())
    {
        FatalErrorIn(args.executable())
            << "No processes for the gases " << gases << " in " << file
            << exit(FatalError);
    }

    Info<< endl;

    const scalar kTe = kB*gasTemperature/eCharge;

    List<swarm> results(fields.size());

    // Method
    const word method(dict.lookupOrDefault<word>("method", "twoTerm"));

    if (method != "twoTerm" && method != "monteCarlo")
    {
        FatalIOErrorIn(args.executable().c_str(), dict)
            << "method must be twoTerm or monteCarlo" << exit(FatalIOError);
    }

    const bool useMonteCarlo = (method == "monteCarlo");

    monteCarloSettings mc;

    mc.nElectrons = dict.lookupOrDefault<label>("nElectrons", 1000);
    mc.gasDensity = dict.lookupOrDefault<scalar>("gasDensity", 3.22e22);
    mc.energySharing = dict.lookupOrDefault<scalar>("energySharing", 0.5);
    mc.equilibrationTimes =
        dict.lookupOrDefault<scalar>("equilibrationTimes", 2);
    mc.sampleTimes = dict.lookupOrDefault<scalar>("sampleTimes", 20);
    mc.minCollisions = dict.lookupOrDefault<scalar>("minCollisions", 2000);
    mc.nBlocks = dict.lookupOrDefault<label>("nBlocks", 10);
    mc.nSpeedCells = dict.lookupOrDefault<label>("nSpeedCells", 4000);

    Random rndGen(dict.lookupOrDefault<label>("seed", 1234));

    List<swarm> twoTermResults(fields.size());
    List<swarm> errors(fields.size());

    mkDir(outputDir);

    scalar epsMax = 1.0;

    forAll(fields, fieldI)
    {
        const scalar EN = fields[fieldI]*Td;

        scalarField F(nCells, 0.0);

        // Adapt the upper end of the grid: the distribution should fall by
        // about ten decades over it
        for (label iter = 0; iter < 40; iter++)
        {
            solveDistribution(processes, gases, EN, kTe, epsMax, F);

            const scalar Fmax = max(mag(F));
            const scalar Fend = mag(F[nCells - 1]);

            if (Fend > 1e-9*Fmax)
            {
                epsMax *= 1.5;
            }
            else
            {
                // Energy at which the distribution has fallen ten decades
                label last = nCells - 1;

                while (last > 0 && mag(F[last]) < 1e-10*Fmax)
                {
                    last--;
                }

                if (last < 0.6*nCells)
                {
                    epsMax *= max((last + 1.0)/(0.8*nCells), 0.3);
                }
                else
                {
                    break;
                }
            }
        }

        const scalar dEps = epsMax/nCells;

        swarm& r = results[fieldI];

        r.EN = EN;
        r.epsMax = epsMax;
        r.meanEnergy = 0;
        r.DN = 0;
        r.DEpsN = 0;
        r.muN = 0;
        r.muEpsN = 0;
        r.rates.setSize(processes.size());
        r.rates = 0.0;

        forAll(F, j)
        {
            const scalar eps = (j + 0.5)*dEps;

            scalar sigmaM, sigmaEps;
            mixtureCrossSections(processes, gases, eps, sigmaM, sigmaEps);

            r.meanEnergy += Foam::pow(eps, 1.5)*F[j]*dEps;
            r.DN += gamma/3*eps/max(sigmaM, 1e-30)*F[j]*dEps;
            r.DEpsN += gamma/3*sqr(eps)/max(sigmaM, 1e-30)*F[j]*dEps;

            forAll(processes, i)
            {
                r.rates[i] += gamma*eps*processes[i](eps)*F[j]*dEps;
            }
        }

        for (label f = 1; f < nCells; f++)
        {
            const scalar eps = f*dEps;

            scalar sigmaM, sigmaEps;
            mixtureCrossSections(processes, gases, eps, sigmaM, sigmaEps);

            r.muN -= gamma/3*eps/max(sigmaM, 1e-30)*(F[f] - F[f - 1]);
            r.muEpsN -= gamma/3*sqr(eps)/max(sigmaM, 1e-30)*(F[f] - F[f - 1]);
        }

        r.DEpsN /= max(r.meanEnergy, 1e-30);
        r.muEpsN /= max(r.meanEnergy, 1e-30);

        Info<< "E/N " << fields[fieldI] << " Td: mean energy "
            << r.meanEnergy << " eV, mu N " << r.muN << " 1/(V m s), D N "
            << r.DN << " 1/(m s), drift velocity " << r.muN*EN
            << " m/s, grid to " << epsMax << " eV" << endl;

        OStringStream name;
        name<< "distribution_" << fields[fieldI] << ".dat";

        OFstream os(outputDir/name.str());

        os  << "# E/N " << fields[fieldI] << " Td: energy [eV]" << tab
            << "F0 [eV^-3/2]" << nl;

        forAll(F, j)
        {
            os  << (j + 0.5)*dEps << tab << F[j] << nl;
        }

        twoTermResults[fieldI] = r;

        if (useMonteCarlo)
        {
            scalarField distribution;
            scalar dV;

            monteCarlo
            (
                processes,
                gases,
                gasTemperature,
                F,
                dEps,
                mc,
                rndGen,
                r,
                errors[fieldI],
                distribution,
                dV
            );

            OStringStream mcName;
            mcName<< "distributionMonteCarlo_" << fields[fieldI] << ".dat";

            OFstream mcOs(outputDir/mcName.str());

            mcOs<< "# E/N " << fields[fieldI] << " Td: energy [eV]" << tab
                << "F0 [eV^-3/2]" << nl;

            forAll(distribution, j)
            {
                mcOs<< sqr((j + 0.5)*dV/gamma) << tab << distribution[j]
                    << nl;
            }
        }
    }


    // Tables against E/N

    if (useMonteCarlo)
    {
        writeTransport
        (
            outputDir/"transport.dat",
            results,
            processes.size(),
            "Monte Carlo; muEps N and DEps N are those of the two-term"
            " solution"
        );

        writeTransport
        (
            outputDir/"transportErrors.dat",
            errors,
            processes.size(),
            "standard errors of the Monte Carlo results"
        );

        writeTransport
        (
            outputDir/"transportTwoTerm.dat",
            twoTermResults,
            processes.size(),
            "two-term approximation"
        );
    }
    else
    {
        writeTransport
        (
            outputDir/"transport.dat",
            results,
            processes.size(),
            "two-term approximation"
        );
    }


    // Tables against the electron temperature, in the formats of the fluid
    // solver. They need Te to increase with E/N: points where it does not
    // are left out

    DynamicList<label> order;
    scalar lastTe = -1;

    forAll(results, k)
    {
        const scalar Te = 2.0/3.0*results[k].meanEnergy*eCharge/kB;

        if (Te > lastTe*(1 + 1e-6))
        {
            order.append(k);
            lastTe = Te;
        }
    }

    const scalar kmol = 6.02214076e26;

    {
        OFstream mu(outputDir/"mu_electron");
        OFstream D(outputDir/"D_electron");

        mu  << "(" << nl;
        D   << "(" << nl;

        forAll(order, j)
        {
            const swarm& r = results[order[j]];
            const scalar Te = 2.0/3.0*r.meanEnergy*eCharge/kB;

            mu  << Te << ' ' << r.muN << nl;
            D   << Te << ' ' << r.DN << nl;
        }

        mu  << ")" << nl;
        D   << ")" << nl;
    }

    forAll(processes, i)
    {
        OStringStream name;
        name<< "rate_" << i;

        OFstream os(outputDir/name.str());

        os  << "(" << nl;

        forAll(order, j)
        {
            const swarm& r = results[order[j]];
            const scalar Te = 2.0/3.0*r.meanEnergy*eCharge/kB;

            os  << "(" << Te << ' '
                << Foam::log(max(r.rates[i]*kmol, 1e-50*kmol)) << ")" << nl;
        }

        os  << ")" << nl;
    }


    // Complete existing tables of the fluid solver

    if (dict.found("extend"))
    {
        const dictionary& ext = dict.subDict("extend");

        // (file, quantity): -1 mobility, -2 diffusion, >= 0 rate of process
        DynamicList<Tuple2<fileName, label> > jobs;

        if (ext.found("mobility"))
        {
            jobs.append
            (
                Tuple2<fileName, label>(fileName(ext.lookup("mobility")), -1)
            );
        }

        if (ext.found("diffusion"))
        {
            jobs.append
            (
                Tuple2<fileName, label>(fileName(ext.lookup("diffusion")), -2)
            );
        }

        if (ext.found("rates"))
        {
            const List<Tuple2<fileName, label> > rates(ext.lookup("rates"));

            forAll(rates, i)
            {
                jobs.append(rates[i]);
            }
        }

        forAll(jobs, jobI)
        {
            fileName table(jobs[jobI].first());
            const label quantity = jobs[jobI].second();

            if (table[0] != '/')
            {
                table = runTime.path()/table;
            }

            // Existing points: all the numbers of the file, in pairs
            DynamicList<scalar> numbers;
            {
                IFstream is(table);

                while (is.good())
                {
                    string line;
                    is.getLine(line);

                    forAll(line, c)
                    {
                        if (line[c] == '(' || line[c] == ')')
                        {
                            line[c] = ' ';
                        }
                    }

                    IStringStream ls(line);

                    while (true)
                    {
                        token t(ls);

                        if (!t.good() || !t.isNumber())
                        {
                            break;
                        }

                        numbers.append(t.number());
                    }
                }
            }

            const label nOld = numbers.size()/2;

            if (nOld < 1)
            {
                WarningIn(args.executable())
                    << "No points in " << table << "; not completed" << endl;
                continue;
            }

            const scalar TeFirst = numbers[0];
            const scalar TeLast = numbers[2*(nOld - 1)];

            const bool rateFormat = (quantity >= 0);

            mv(table, table + ".orig");

            OFstream os(table);

            os  << "(" << nl;

            label nBelow = 0;
            label nAbove = 0;

            // Computed points below, the existing points, computed above
            for (label pass = 0; pass < 3; pass++)
            {
                if (pass == 1)
                {
                    for (label j = 0; j < nOld; j++)
                    {
                        if (rateFormat)
                        {
                            os  << "(" << numbers[2*j] << ' '
                                << numbers[2*j + 1] << ")" << nl;
                        }
                        else
                        {
                            os  << numbers[2*j] << ' ' << numbers[2*j + 1]
                                << nl;
                        }
                    }

                    continue;
                }

                forAll(order, j)
                {
                    const swarm& r = results[order[j]];
                    const scalar Te = 2.0/3.0*r.meanEnergy*eCharge/kB;

                    if
                    (
                        (pass == 0 && Te < TeFirst*(1 - 1e-6))
                     || (pass == 2 && Te > TeLast*(1 + 1e-6))
                    )
                    {
                        const scalar value =
                        (
                            quantity == -1 ? r.muN
                          : quantity == -2 ? r.DN
                          : Foam::log(max(r.rates[quantity]*kmol, 1e-50*kmol))
                        );

                        if (rateFormat)
                        {
                            os  << "(" << Te << ' ' << value << ")" << nl;
                        }
                        else
                        {
                            os  << Te << ' ' << value << nl;
                        }

                        (pass == 0 ? nBelow : nAbove)++;
                    }
                }
            }

            os  << ")" << nl;

            Info<< "Completed " << table << ": " << nOld
                << " points kept (" << TeFirst << " to " << TeLast
                << " K), " << nBelow << " added below and " << nAbove
                << " above; original saved as .orig" << endl;
        }
    }

    Info<< nl << "Tables written to " << outputDir << nl << nl << "End" << nl
        << endl;

    return 0;
}


// ************************************************************************* //
