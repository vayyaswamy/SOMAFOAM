	ionFlux = mspm().totalIonFlux();

	mspm().correct(Y, E, chemistry, fields);

	pem.ecorrect(chemistry, E);

	pem.icorrect(chemistry);

	pem.gcorrect(chemistry, E);

	thermo.correct();

	rho = thermo.rho();

	volVectorField F = mspm().netChargeFlux();

	if (emcModelRef == "power")
	{
		ddtE = fvc::ddt(E);
	}
