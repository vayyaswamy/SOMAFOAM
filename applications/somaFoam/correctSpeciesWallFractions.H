// Mass fractions on fixed-value walls from the densities (start-up).
//
// The species files hold number densities with their boundary conditions
// (they are written from N, so that a restart reads the boundary settings
// back). The mass fractions Y are read from the same files and their cell
// values are converted from N (Y = N W/(rho A)) by the plasma model, but an
// assignment does not change a fixed-value boundary, so such a boundary kept
// the density of the file as a "mass fraction": 0 at a fresh start, the
// density floor (e.g. 1e4) after a restart, which made the mixture of the
// wall faces meaningless and the thermophysical update fail. Here every
// fixed-value boundary of Y is set from the boundary density once at
// start-up; afterwards it stays consistent.
{
    const volScalarField rhoStart(thermo.rho());

    forAll(Y, i)
    {
        const volScalarField& Ni = mspm().N(i);

        forAll(Y[i].boundaryField(), patchi)
        {
            fvPatchScalarField& Yp = Y[i].boundaryField()[patchi];

            if (Yp.fixesValue() && Yp.size())
            {
                Yp == Ni.boundaryField()[patchi]*mspm().W(i)
                   /rhoStart.boundaryField()[patchi]/plasmaConstants::A.value();
            }
        }
    }
}
