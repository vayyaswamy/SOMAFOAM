
	
	ionFlux = mspm().totalIonFlux();
	mspm().correct(Y, E, chemistry, fields);

	


	kappa_e = mspm().electronConductivity(chemistry);

		

	//Info << "Done with electron temperature" << endl;

	thermo.correct();

	rho = thermo.rho();

	gradTe = mspm().gradTe();



	pem.icorrect(chemistry);

	//Info << "About to gcorrect " << endl;

	pem.gcorrect(chemistry, E);

	

	Jnet = mspm().netChargeFlux();

	// Compute total current density (conduction + displacement) in plasma region
	Jtot = Jnet + eps*(E - E.oldTime())/runTime.deltaT().value();

	if (emcModelRef == "power")
	{
		ddtE = fvc::ddt(E);
	}
