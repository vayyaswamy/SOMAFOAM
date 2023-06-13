	
	ionFlux = mspm().totalIonFlux();
	mspm().correct(Y, E, chemistry, fields);

	//Info << "Densities corrected " << endl;

	kappa_e = mspm().electronConductivity(chemistry);

	//Info << "kappe_e done " << endl;

	pem.ecorrect(chemistry, E);	

	//Info << "Te corrected " << endl;

	thermo.correct();


	gradTe = mspm().gradTe();

	

	//Info << "About to icorrect " << endl;

	pem.icorrect(chemistry);

	//Info << "About to gcorrect " << endl;

	pem.gcorrect(chemistry, E);

	thermo.correct();

	rho = thermo.rho();

	Jnet = mspm().netChargeFlux();

	// Compute total current density (conduction + displacement) in plasma region
	Jtot = Jnet + eps*(E - E.oldTime())/runTime.deltaT().value();

	if (emcModelRef == "power")
	{
		ddtE = fvc::ddt(E);
	}
