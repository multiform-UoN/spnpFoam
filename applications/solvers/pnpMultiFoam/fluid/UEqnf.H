//Momentum predictor
tmp<fvVectorMatrix> tUEqn
(
    //Stokes flow approximation
    - fvm::laplacian(mu,U)
    + fvc::reconstruct
      (
         fvc::interpolate(rho)*fvc::snGrad(V)*mesh.magSf()
      )
    ==
    unitDensity*fvOptions(U)
);

fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();
fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    //Solve momentum predictor
    solve
    (
        UEqn + fvc::reconstruct(fvc::snGrad(p)*mesh.magSf())
    );

    fvOptions.correct(U);
}
