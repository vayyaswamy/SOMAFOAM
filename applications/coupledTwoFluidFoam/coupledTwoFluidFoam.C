/*----------------------:-----------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    coupledTwoFluidFoam

Description
    Solver for basic plasma fluid model with electron and one ion species fully-coupled
    with Poisson's equation.
    There is no electron energy equation and no collisions included.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvBlockMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting iteration loop\n" << endl;

    //scalar e = 1.602E-19;
    //scalar epsilon0 = 8.854E-12;
    //
    dimensionedScalar e 
    (
        "e",
        dimensionSet(0,0,1,0,0,1,0),
        scalar(1.602e-19)
    );

    dimensionedScalar epsilon0 
    (
        "epsilon0",
        dimensionSet(-1,-3,4,0,0,2,0),
        scalar(8.854e-12)
    );

    




    while (runTime.loop())
    {
        Info<< "Iteration = " << runTime.timeName() << nl << endl;

        // Initialize the plasma block system (matrix, source and reference to Up)
        fvBlockMatrix<vector9> plasmaEqn(plasma);

	   // Continuity equation for electrons
	   fvScalarMatrix neEqn
	   (
	       fvm::ddt(ne) + fvm::div(neFlux,ne) - Se 
	   );

	   // Momentum equation for electrons
        fvVectorMatrix UeEqn
	   (
	       fvm::ddt(ne,Ue) + fvm::div(UeFlux,Ue) + e*ne*E/me + e*fvc::grad(ne*Te)/me 
	   );


	   // Continuity equation for ion 1
        fvScalarMatrix n1Eqn
	   (
 	      fvm::ddt(n1) + fvm::div(n1Flux,n1) - S1 
	   );

	   // Momentum equation for ion 1
	   fvVectorMatrix U1Eqn
	   (
 	      fvm::ddt(n1,U1) + fvm::div(U1Flux,U1) - e*n1*E/m1 + e*fvc::grad(n1*T1)/m1 
	   );

        // Poisson's equation
	   fvScalarMatrix phiEqn
	   (
	    fvm::laplacian(phi) + e*(n1 - ne)/epsilon0 
	   );


       


       // Creating auxillary fields from solved variables
        Gammae = ne*Ue;
	    Gamma1 = n1*U1;

        // Update electron flux
	   neFlux = fvc::interpolate(Ue) & mesh.Sf();
	   UeFlux = fvc::interpolate(ne)*neFlux;

	   // Update ion 1 flux
	   n1Flux = fvc::interpolate(U1) & mesh.Sf();
	   U1Flux = fvc::interpolate(n1)*n1Flux;
       
        // Calculate E-field
        E = -fvc::grad(phi);


        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
