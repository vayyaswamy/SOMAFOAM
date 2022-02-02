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
    plasmaFluidFoam

Description
    Solver for basic plasma fluid model with electron and two ion species.
    There is no electron energy equation and no collisions included.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

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

	// Continuity equation for electrons
	solve
	(
	    fvm::ddt(ne) + fvm::div(neFlux,ne) - Se
	);

	// Momentum equation for electrons
        solve
	(
	    fvm::ddt(ne,Ue) + fvm::div(UeFlux,Ue) + e*ne*E/me + e*fvc::grad(ne*Te)/me
	);


	// Continuity equation for ion 1
        solve
	(
 	    fvm::ddt(n1) + fvm::div(n1Flux,n1) - S1 
	);

	// Momentum equation for ion 1
	solve
	(
 	    fvm::ddt(n1,U1) + fvm::div(U1Flux,U1) - e*n1*E/m1 + e*fvc::grad(n1*T1)/m1 
	);


	// Continuity equation for ion 2
	solve
	(
            fvm::ddt(n2) + fvm::div(n2Flux,n2) - S2
	);

	// Momentum equation for ion 2
	solve
	(
	    fvm::ddt(n2,U2) + fvm::div(U2Flux,U2) - e*n2*E/m2 +e*fvc::grad(n2*T2)/m2
	);

        // Poisson's equation
	solve
	(
	    fvm::laplacian(phi) + e*(n1 + n2 - ne)/epsilon0
	);

        Gammae = ne*Ue;

	Gamma1 = n1*U1;

	Gamma2 = n2*U2;

        // Update electron flux
	neFlux = fvc::interpolate(Ue) & mesh.Sf();
	UeFlux = fvc::interpolate(ne)*neFlux;

	// Update ion 1 flux
	n1Flux = fvc::interpolate(U1) & mesh.Sf();
	U1Flux = fvc::interpolate(n1)*n1Flux;
       
        // Update ion 2 flux
	n2Flux = fvc::interpolate(U2) & mesh.Sf();
	U2Flux = fvc::interpolate(n2)*n2Flux;

        // Calculate E-field
        E = -fvc::grad(phi);

	ni = n1 + n2;

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
