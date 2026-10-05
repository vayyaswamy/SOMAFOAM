

    ionFlux = mspm().totalIonFlux();
    mspm().correct(Y, E, chemistry, fields);


    kappa_e = mspm().electronConductivity(chemistry);


    thermo.correct();

    rho = thermo.rho();

    gradTe = mspm().gradTe();


    pem.icorrect(chemistry);
    mspm().limitField(thermo.Tion());

    Info << "About to gcorrect " << endl;

    pem.gcorrect(chemistry, E);
    mspm().limitField(thermo.T());


    Jnet = mspm().netChargeFlux();

    // Compute total current density (conduction + displacement) in plasma region
    Jtot = Jnet + eps*(E - E.oldTime())/runTime.deltaT().value();

    if (emcModelRef == "power")
    {
        ddtE = fvc::ddt(E);
    }
