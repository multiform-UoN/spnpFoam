// Set current mesh to fluid region i
fvMesh& mesh(fluidRegions[i]);

//Set current fields to fluid region i
volScalarField& V(VFluid[i]);
volScalarField& T(TFluid[i]);
surfaceScalarField& phiNP(phiNPFluid[i]);
volVectorField& U(UFluid[i]);
volScalarField& rho(rhoFluid[i]);
volScalarField& p(pFluid[i]);
surfaceScalarField& phi(phiFluid[i]);

// Set current concentrations to fluid region i
PtrList<volScalarField>& cPtrL(cPtrLFluid[i]);

// Set species dictionary to fluid region i
PtrList<dictionary>& species(speciesFluid[i]);

// Set effective dielectric to fluid region i
dimensionedScalar& epsilon_(epsilon_Fluid[i]);

// Set dynamic viscocity to fluid region i
dimensionedScalar& mu(mu_Fluid[i]);

// Set electro-neutrality and diffusion dominant flags to fluid region i
const bool&	electroNeutrality(electroNeutralityFluid[i]);

// Set reference cell and value to fluid region i
const label& VRefCell(VRefCellFluid[i]);
const scalar& VRefValue(VRefValueFluid[i]);
const label& pRefCell(pRefCellFluid[i]);
const scalar& pRefValue(pRefValueFluid[i]);

// // Set fvOptions to fluid region i
fv::options& fvOptions(fvOptionsFluid[i]);

// Set list of Max Co #'s for each species to fluid region i
List<scalar>& spCoNum(speciesCoNumFluid[i]);

// Set list of Mean Co #'s for each species to fluid region i
List<scalar>& spMeanCoNum(speciesMeanCoNumFluid[i]);

// Find # of loop controls set for fluid region i
pimpleNoLoopControl& pimple(pimples.pimple(i));
