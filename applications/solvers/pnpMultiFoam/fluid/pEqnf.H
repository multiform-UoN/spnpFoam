// Store reciprocal of central coeffs of UEqn fvMatrix
volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));


// Store off-diagnoal components of UEqn fvMatrix, times by reciprocal central coeffs
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p)); // U = H/A - presgrad/A - bodyforce/A

// Determine phiHbyA = sumOverFaces((HbyA)_f . S_f)
surfaceScalarField phiHbyA
(
    "phiHbyA",
    fvc::flux(HbyA)
);

if (p.needReference()) // CHECK IF phiV needed
{
    fvc::makeRelative(phiHbyA, U);
    adjustPhi(phiHbyA, U, p);
    fvc::makeAbsolute(phiHbyA, U);
}

if (pimple.nCorrPiso() <= 1)
{
    tUEqn.clear();
}

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, U, phiHbyA, rAUf);

// Non-orthogonal pressure corrector loop
// Pressure correction eqn: div(grad(p)/A) = div((H/A) - (rho*grad(V)/A))
while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pEqn
    (
         - fvm::laplacian(rAU, p)
        ==
         - fvc::div(phiHbyA)
    );

    pEqn.setReference(pRefCell, pRefValue);
    pEqn.relax();
    pEqn.solve();

    // Explicitly relax pressure for momentum corrector
    p.relax();

    if (pimple.finalNonOrthogonalIter())
    {
        phi = phiHbyA + pEqn.flux();
        U = HbyA + rAU*fvc::reconstruct(pEqn.flux()/rAUf);
        U.correctBoundaryConditions();
        fvOptions.correct(U);
    }
}

#include "continuityErrs.H"
