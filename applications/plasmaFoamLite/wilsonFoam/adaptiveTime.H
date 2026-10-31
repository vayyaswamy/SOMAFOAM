    surfaceScalarField magPhi = mag(fvc::interpolate(EField*mobility[0]) & mesh.Sf());
    if (transportModel[0]=="driftDiffusion")
        surfaceScalarField magPhi = mag(fvc::interpolate(EField*mobility[0]) & mesh.Sf());
    else
        surfaceScalarField magPhi = mag(fvc::interpolate(velocity[0]) & mesh.Sf());
    double plasmaFreq = Foam::sqrt(gMax(nDensity[0])*1.6e-19*1.6e-19/9.1e-31/8.85e-12)*6.238;//cyclespersecond
    surfaceScalarField SfUfbyDelta =
        mesh.surfaceInterpolation::deltaCoeffs()*magPhi;
    if (average(SfUfbyDelta/mesh.magSf()).value()!=0){
        double dt=maxCo/(average(SfUfbyDelta/mesh.magSf()).value());
        //dt=std::min(dt,1/plasmaFreq);
        dt=std::min(dt,deltaTMax);
        dt=std::max(dt,deltaTMin);
        Info << "predicted by plasma frequency"<< 1/plasmaFreq << endl;
        Info << "max value of electronDensity" << gMax(nDensity[0]) << endl;
        runTime.setDeltaT(dt*dimensionCorrectionTime);
    }
    Info << "Delta Time::::= "<<runTime.deltaT().value() << endl;