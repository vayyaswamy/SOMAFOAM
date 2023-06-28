	int icorr = 1;

	scalar initialResidual = 1.0;

	while ((icorr++ <= 5) && (initialResidual >= 1e-5))
	{
		ionFlux = mspm().totalIonFlux();
		initialResidual = mspm().correct(Y, E, chemistry, fields);

	}

	//Info << "Done with number density" << endl;
	

	//ionFlux = mspm().totalIonFlux();
	//mspm().correct(Y, E, chemistry, fields);	
	

	kappa_e = mspm().electronConductivity(chemistry);

	pem.ecorrect(chemistry, E);	

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
