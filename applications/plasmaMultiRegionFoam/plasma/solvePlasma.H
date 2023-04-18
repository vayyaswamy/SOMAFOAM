	ionFlux = mspm().totalIonFlux();

	mspm().correct(Y, E, chemistry, fields);

	kappa_e = mspm().electronConductivity(chemistry);

	pem.ecorrect(chemistry, E); 

	pem.icorrect(chemistry); 

	pem.gcorrect(chemistry, E);

	gradTe = mspm().gradTe();

	thermo.correct();

	rho = thermo.rho();

	volVectorField F = mspm().netChargeFlux();

	conductionCurrent = mspm().netChargeFlux();

	totalCurrent = mspm().netChargeFlux() + plasmaConstants::epsilon0*fvc::ddt(E);

	if (emcModelRef == "power")
	{
		ddtE = fvc::ddt(E);
	}

	{
		forAll(dielectricRegions, i)
		{
			label patchID = mesh.boundaryMesh().findPatchID(rp.dielectricRegionNames()[i]);			
			
			// Euler scheme hardcoded
			scalarField FN = F.boundaryField()[patchID].patchInternalField() & mesh.Sf().boundaryField()[patchID] / mesh.magSf().boundaryField()[patchID];

			scalarField FSN = F.boundaryField()[patchID].snGrad() & mesh.Sf().boundaryField()[patchID] / mesh.magSf().boundaryField()[patchID];

			surfC.boundaryField()[patchID] =
				surfC.oldTime().boundaryField()[patchID]
				+ eCharge.value()*(FN + (FSN / mesh.boundary()[patchID].deltaCoeffs()))*runTime.deltaT().value();
		}
	}
