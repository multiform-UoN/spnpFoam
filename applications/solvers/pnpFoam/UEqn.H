MRF.correctBoundaryVelocity(U);  //added

// Stokes equation
tmp<fvVectorMatrix> tUEqn
(
    -
    fvm::laplacian(mu,U)
    +
    bodyForce*rho*fvc::grad(V)
    ==
    unitDensity*fvOptions(U)
);

fvVectorMatrix& UEqn = tUEqn.ref();
  
UEqn.relax();
  
fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    //Solve momentum predictor
    solve(UEqn == -fvc::grad(p));

    fvOptions.correct(U);
}