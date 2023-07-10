/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    emFoam

Description
    Solver for (quasi-)stationary electromagnetic fields produced by steady or 
    alternating currents with a given frequency. The vector potential is 
    solved in the frequency domain using a block coupled matrix solver so that 
    even high-frequency electromagnetic phenomena can be treated well. 

    Note: "Jcoil" refers to the source current producing the electromagnetic 
    field I kept this naming convention to be consistent with the documentation 
    in my master thesis, where you can read more about the theory behind the 
    solver: 
    Busse, Christian: Numerical Modeling of an Inductively Coupled 
    Plasma (ICP). Ilmenau 2019.  https://doi.org/10.22032/dbt.40314

    In the case setup the following inputs need to be specified
    in "../case/0":
        - Source current density Jcoil in A/m^2
        - Electrical conductivity sigma in A/Vm
    in "../case/constant/physicalProperties":
        - Magnetic permeability muMag (default muMag=mu0) in Vs/Am
        - Current frequency w in 1/s

    The output quantities are:
        - Magnetic vector potential A in Vs/m
        - Magnetic flux density B in Vs/m^2
        - Magnetic field strength H in A/m
        - Induced current density Jind in A/m^2
        - (Time-averaged) Joule heat density in W/m^3  
        - (Time-averaged) Lorentz-force density fL in N/m^3

Author
    Christian Busse

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
//#include "coupledFvMatrices.H"
//#include "fvBlockMatrix.H" // load block coupled matrix solver

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



void normalize
    (
    volVectorField& AR,
    volVectorField& AI,
    volVectorField& J,
    const scalar factor
    )
{
    AR *= factor;
    AI *= factor;
    J  *= factor;
}


void report
    (
    const volVectorField& AR,
    const volVectorField& AI,
    const volVectorField& AR_last,   // CQ: pass AR_last to calculate ar_fracMax
    const volVectorField& AI_last,   // CQ: pass AI_last to calculate ai_fracMax
    const volVectorField& J,
    const volScalarField& sigma,
    const dimensionedScalar& w
    )
{
    Info<< "----- Reporting: " << nl << endl;

    const volScalarField Amag = sqrt((AR & AR) + (AI & AI)); // & - dot product of two vectors
    const volScalarField qJ = sigma/2.0 * sqr(w) * sqr(Amag); // Dissipated power density in [W/m^3]   
    dimensionedScalar QJ = fvc::domainIntegrate(qJ);
    Info<< "QJ     : " << QJ.value() << nl
        << "(Compare to value caluclated with emFoam)" << nl << endl;

    // CQ: extra tests
    scalar ar_min = 1e20;
    scalar ar_avg = 0.0;
    scalar ar_max = -1e20;
    scalar ar_fracMax = 0.;   // CQ: maximum fractional change of AR 
    scalar ai_min = 1e20;
    scalar ai_avg = 0.0;
    scalar ai_max = -1e20;
    scalar ai_fracMax = 0.;   // CQ: maximum fractional change of AI
    scalar count = 0.0;
    forAll (AR, cellI)
    {
        ar_min = min(AR.internalField()[cellI][1], ar_min);
        ar_max = max(AR.internalField()[cellI][1], ar_max);
        ar_avg += AR.internalField()[cellI][1];
        ar_fracMax = max(AR.internalField()[cellI][1]/(AR_last.internalField()[cellI][1]+1.e-40), ar_fracMax); // CQ: compute maximum fractional change of AR
        ai_min = min(AI.internalField()[cellI][1], ai_min);
        ai_max = max(AI.internalField()[cellI][1], ai_max);
        ai_avg += AI.internalField()[cellI][1];
        ai_fracMax = max(AI.internalField()[cellI][1]/(AI_last.internalField()[cellI][1]+1.e-40), ai_fracMax); // CQ: compute maximum fractional change of AI
        count += 1.0;
    }
    ar_avg /= count;
    ai_avg /= count;

    Info<< "max(AR): " << ar_max   << nl
        << "min(AR): " << ar_min   << nl
        << "avg(AR): " << ar_avg   << nl
        << "ar_fracMax: " << ar_fracMax << nl        // CQ: export maximum fractional change of AR
        << nl
        << "max(AI): " << ai_max   << nl
        << "min(AI): " << ai_min   << nl
        << "avg(AI): " << ai_avg   << nl
        << "ai_fracMax: " << ai_fracMax << nl        // CQ: export maximum fractional change of AI
        << endl;
}


int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"
#   include "initConvergenceCheck.H"

    bool isFirstPass = true;

    Info<< "\nStarting time loop\n" << endl;
    while (runTime.loop())
    {
/*
#       include "readBlockSolverControls.H" // I didn't use these, but someone might want to have it
#       include "readFieldBounds.H"
*/   
        AR.storePrevIter();   // CQ: store previter so we can relax on AR
        AI.storePrevIter();   // CQ: store previter so we can relax on AI
        Info<< "Time = " << runTime.timeName() << nl << endl;

        DA = w * muMag * sigma;

// CQ: comment "guess" part, don't use both blockCoupled and segregated solver together
//        // 2) Calculate guess
//        //    Use the original block coupled code to calculate a "guess" for
//        //    the solution. To get an ineact guess, the maximum number of
//        //    iterations used to solve the A equation (sest in
//        //    system/fvSolution) should be set to low value, say 50.
//        if (isFirstPass)
//        {
//            Info<< "----- Calculating initial 'guess': " << nl << endl;
//
//            // Initialize the vector potential matrix
//            Info<< "Constructing block coupled matrix" << endl;
//            fvBlockMatrix<vector6> AEqn(A);
//
//                fvVectorMatrix AREqn
//                (
//                    fvm::laplacian(AR)
//                    ==
//                    - muMag * Jcoil
//                ); 
//
//                fvVectorMatrix AIEqn
//                (
//                    fvm::laplacian(AI)
//                );
//
//            //insert fvVectorMatrix equations into the fvBlockMatrix
//            AEqn.insertEquation(0, AREqn);
//            AEqn.insertEquation(3, AIEqn);
//
//            // Add off-diagonal coupling terms
//            AEqn.insertEquationCoupling(0, 3, DA);
//            AEqn.insertEquationCoupling(1, 4, DA);
//            AEqn.insertEquationCoupling(2, 5, DA);
//            AEqn.insertEquationCoupling(3, 0, -DA);
//            AEqn.insertEquationCoupling(4, 1, -DA);
//            AEqn.insertEquationCoupling(5, 2, -DA);
//
//            // Solve the block matrix
//            Info<< "Solving block coupled matrix" << endl;
//            maxResidual = cmptMax(AEqn.solve().initialResidual());
//
//            // Retrieve solution of the vector potential A
//            AEqn.retrieveSolution(0, AR.internalField());
//            AEqn.retrieveSolution(3, AI.internalField());
//
//            report(AR, AI, Jcoil, sigma, w);      
//            AR_last = AR;  // CQ: save guess AR as the first "last" value for AR used in relaxation
//            AI_last = AI;  // CQ: save guess AI as the first "last" value for AI used in relaxation
//
//            Info<< "First guess calculated" << endl;
//            isFirstPass = false;
//        }
//        else                                        // CQ: else if not the first step (and use segregated)
        // 3) Solve for AR and AI using SOR
        //    Note, relaxation parameters are set in system/fvSolution
//        {
            Info<< "----- Solving for this iteration: " << nl << endl;
            Info<< "Solving for AR" << endl;

            volVectorField ARSource = muMag*Jcoil + DA*AI ;
            fvVectorMatrix AREqn               
            ( 
                  fvm::laplacian(AR)

                  + fvm::SuSp((ARSource/AR),AR)
                ==
                0
            );
            //AREqn.relax();
            
            AREqn.solve();

            AR.relax();

            Info<< "Solving for AI" << endl;
            volVectorField AISource = -DA*AR;
            fvVectorMatrix AIEqn
            (
                  fvm::laplacian(AI)
                  +fvm::SuSp((AISource/AI),AI)
                ==
                  0
            );
            //AIEqn.relax();
            
            AIEqn.solve();

            AI.relax();
//        }                                    // CQ: no blockCoupled "if else"

        report(AR, AI, AR_last, AI_last, Jcoil, sigma, w);

        deltaAR = AR - AR_last ;
        deltaAI = AI - AI_last ;

        AR_last = AR;     // CQ: update AR_last to new value  
        AI_last = AI;     // CQ: update AI_last to new value


        #include "convergenceCheck.H"

        Info<< "EM Solver: ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}
