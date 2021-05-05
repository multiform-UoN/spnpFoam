//Momentum predictor
tmp<fvVectorMatrix> tUEqn
(
    //Stokes flow approximation
    -
    fvm::laplacian(mu,U)
    +
    rho*fvc::grad(V)
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
