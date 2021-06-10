{
  //-Update charge density
  rho *= scalar(0); //- Simply reset to zero

  if (!electroNeutrality){

    forAll(species,sp)
    {
      dimensionedScalar Z(species[sp].lookup("Z"));
      rho += Fy_*cPtrL[sp]*Z;
    }

  }

  // If flow isn't frozen
  if(!pimple.frozenFlow())
  {
    // Solve Stokes
    if(faceMomentum)
    {
        #include "UEqnf.H"

        // Solve pressure correction
        while(pimple.correct())
        {
          #include "pEqnf.H"
        }

    }
    else
    {
        #include "UEqn.H"

        // Solve pressure correction
        while(pimple.correct())
        {
          #include "pEqn.H"
        }
    }
  }

  //- Solve Poisson
  while (pimple.correctNonOrthogonal())
  {
    fvScalarMatrix VEqn
    (
      - fvm::laplacian(epsilon_,V)
      ==
      rho
    );

    // Set the reference value for V to RefCell & RefValue, do nothing if not
    VEqn.setReference(VRefCell, VRefValue, forceRef);

    VEqn.relax();
    VEqn.solve();
  }


  //- Solve species
  forAll(species,sp)
  {

    //- Collect pointers and coefficients
    volScalarField&   C(cPtrL[sp]);
    surfaceScalarField& phiNP(phiPtrL[sp]);
    dimensionedScalar D(species[sp].lookup("D"));
    dimensionedScalar Z(species[sp].lookup("Z"));

    //- Update Nerst-Planck flux to specific concentration
    phiNP = phi - e_*fvc::snGrad(V)*mesh.magSf()/(k_*fvc::interpolate(T))*D*Z;

    //Check if electro-neutrality is wanted
    if (electroNeutrality && sp == species.size()-1)
    {
      //Look up the first species
      const dimensionedScalar& Z_ini(species[0].lookup("Z"));
      const volScalarField& C_ini(cPtrL[0]);

      //Set concentration C to first species
      C  = -Z_ini*C_ini/Z;

      //Loop through all species except final species
      for (int i=1;i <= species.size()-2;i++)
      {
        //Look up the next species
        const dimensionedScalar& Z_iter(species[i].lookup("Z"));
        const volScalarField& C_iter(cPtrL[i]);

        //Add next species to C
        C -=  Z_iter*C_iter/Z;

      }

      C.correctBoundaryConditions();

    }
    else
    {

      //- Non-orthogonal correction loop
      while (pimple.correctNonOrthogonal())
      {

        C.storePrevIter();
        fvScalarMatrix CEqn
        (
          fvm::ddt(C)
          + fvm::div(phiNP,C,"div(phiNP,C)")
          - fvm::laplacian(D,C,"laplacian(D,C)")
        );

        //constrainFluxes(CEqn); // only for codedFlux

        CEqn.relax();
        CEqn.solve();

      }
    }

    // Update tot mass
    speciesTotConcFluid[i][sp] =
        Foam::gSum(cPtrLFluid[i][sp]().field()*fluidRegions[i].V());

  }
}
