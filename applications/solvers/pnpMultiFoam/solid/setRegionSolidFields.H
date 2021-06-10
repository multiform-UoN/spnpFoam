fvMesh& mesh = solidRegions[i];

//Set current fields to solid region i
volScalarField& V = VSolid[i];
volScalarField& T = TSolid[i];
volScalarField& rho = rhoSolid[i];

//Set current concentrations to solid region i
PtrList<volScalarField>& cPtrL = cPtrLSolid[i];

//Set species dictionary to solid region i
PtrList<dictionary>& species = speciesSolid[i];

//Set effective dielectric to solid region i
dimensionedScalar& epsilon_ = epsilon_Solid[i];

//Set electro-neutrality and diffusion dominant flags to fluid region i
const bool&	electroNeutrality = electroNeutralitySolid[i];

//Set reference cell and value to fluid region i
const label& VRefCell = VRefCellSolid[i];
const scalar& VRefValue = VRefValueSolid[i];
const bool forceRef = solidRegions[i].solutionDict().subDict("PIMPLE").lookupOrDefault<bool>("forceVRef",false);

// Set list of Max Co #'s for each species to fluid region i
List<scalar>& spDiffNum(speciesDiffNumSolid[i]);

// Set list of Mean Co #'s for each species to fluid region i
List<scalar>& spMeanDiffNum(speciesMeanDiffNumSolid[i]);

//Find # of loop controls set for solid region i
solidNoLoopControl& pimple = pimples.solid(i);
