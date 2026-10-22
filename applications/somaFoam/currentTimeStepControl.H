// Smallest time step suggested by the electrodes (plasmaPotential patches)
// that have the current-based timeStepControl option on; GREAT if none.
namespace Foam
{

scalar currentLimitedDeltaT(const volScalarField& Phi)
{
    scalar dt = GREAT;

    forAll(Phi.boundaryField(), patchi)
    {
        const fvPatchScalarField& pf = Phi.boundaryField()[patchi];

        if (isA<plasmaPotential>(pf))
        {
            const plasmaPotential& pp = refCast<const plasmaPotential>(pf);

            if (pp.timeStepControl())
            {
                dt = min(dt, pp.suggestedDeltaT());
            }
        }
    }

    return returnReduce(dt, minOp<scalar>());
}

}
