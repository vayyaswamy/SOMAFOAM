	
	while (pimple.correct()){
	ionFlux = mspm().totalIonFlux();
	mspm().correct(Y, E, chemistry, fields);
	}
	

	kappa_e = mspm().electronConductivity(chemistry);

	while (pimple.correct()){
	pem.ecorrect(chemistry, E);
	}
		
	thermo.correct();
	gradTe = mspm().gradTe();

	pem.icorrect(chemistry);
	
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
