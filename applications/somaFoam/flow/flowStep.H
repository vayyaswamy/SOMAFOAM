// Gas flow (createFlowFields.H), after each plasma step: a transient flow
// takes one rhoPimpleFoam step with the plasma force and heating of this
// plasma step; a steady flow with updateInterval is solved again with their
// averages over that interval. A frozen flow, and a steady flow between
// updates, are held.
if (flowTransient || (flowSteady && flowUpdateInterval > 0))
{
    // The mesh changes with adaptive refinement: the sums start again
    if (plasmaForceSum().size() != mesh.nCells())
    {
        plasmaForceSum().setSize(mesh.nCells());
        plasmaHeatSum().setSize(mesh.nCells());
        plasmaForceSum() = vector::zero;
        plasmaHeatSum() = 0;
        flowTimeSum = 0;
    }

    {
        const scalar dtPlasma = runTime.deltaTValue();

        // Force of the field on the space charge, e sum_i z_i N_i E
        const scalarField charge(mspm().netCharge()().internalField());

        plasmaForceSum() +=
            dtPlasma*plasmaConstants::eCharge.value()*charge
           *E.internalField();

        if (flowEnergy)
        {
            plasmaHeatSum() +=
                dtPlasma*mspm().ionTempSource(chemistry, E)().internalField();
        }

        flowTimeSum += dtPlasma;
    }

    if
    (
        flowTransient
     || flowTimeSum >= flowUpdateInterval*(1 - SMALL)
    )
    {
        plasmaForcePtr().internalField() = plasmaForceSum()/flowTimeSum;
        plasmaForcePtr().correctBoundaryConditions();
        plasmaHeatPtr().internalField() = plasmaHeatSum()/flowTimeSum;
        plasmaHeatPtr().correctBoundaryConditions();

        if (flowTransient)
        {
            #include "flowTransient.H"
        }
        else
        {
            Info<< "Steady flow with the plasma force and heating averaged"
                << " over " << flowTimeSum << " s" << endl;

            #include "flowSteady.H"
        }

        plasmaForceSum() = vector::zero;
        plasmaHeatSum() = 0;
        flowTimeSum = 0;
    }
}
