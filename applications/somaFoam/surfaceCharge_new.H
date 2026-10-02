 {
    forAll(Phi.boundaryField(), patch)
    {

     if (Phi.boundaryField()[patch].type() == "dielectricSideWall")
     {
       surfC.boundaryField()[patch] = surfC.oldTime().boundaryField()[patch] + (Jnet.boundaryField()[patch] & mesh.Sf().boundaryField()[patch]) / mesh.magSf().boundaryField()[patch] *runTime.deltaT().value();
     }
     if (Phi.boundaryField()[patch].type() == "externalCircuitPotential")
     {
         scalar Q = gSum(Jnet.boundaryField()[patch] & mesh.Sf().boundaryField()[patch]);
         capCharge.boundaryField()[patch] = capCharge.oldTime().boundaryField()[patch] + Q*runTime.deltaT().value() ;
    }
}
}
