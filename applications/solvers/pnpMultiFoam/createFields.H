
Info<< "Reading electrokineticProperty files\n" << endl;

//Pointer lists for electroKineticProperty dictionaries
PtrList<IOdictionary> SolidElectroKineticProperties(solidRegions.size());
PtrList<IOdictionary> FluidElectroKineticProperties(fluidRegions.size());

//Pointer lists for relative and effective dielectrics and dynamic viscocity per region
PtrList<dimensionedScalar> epsilonR_Solid(solidRegions.size());
PtrList<dimensionedScalar> epsilon_Solid(solidRegions.size());

PtrList<dimensionedScalar> epsilonR_Fluid(fluidRegions.size());
PtrList<dimensionedScalar> epsilon_Fluid(fluidRegions.size());
PtrList<dimensionedScalar> mu_Fluid(fluidRegions.size());


// Run through all solid regions
forAll(solidRegions, i){

	// Set ith electrokineticProperties to dictionary of file
	SolidElectroKineticProperties.set(

		i,
		new IOdictionary(

			IOobject(

				"electrokineticProperties",
				runTime.constant(),
				solidRegions[i],
				IOobject::MUST_READ_IF_MODIFIED,
				IOobject::NO_WRITE
			)
		)
	);

// Construct dielectric values from dictionaries
//(CHANGED TO NON-CONSTANT TO ALLOW VARYING VALUES BASED ON CONCENTRATION,ETC!!!!)
epsilonR_Solid.set(i,new dimensionedScalar(SolidElectroKineticProperties[i].lookup("epsilonR")));
epsilon_Solid.set(i,new dimensionedScalar(epsilonR_Solid[i] * epsilonZero_));

}

forAll(fluidRegions, i){

	FluidElectroKineticProperties.set(

		i,
		new IOdictionary(

			IOobject(

				"electrokineticProperties",
				runTime.constant(),
				fluidRegions[i],
				IOobject::MUST_READ_IF_MODIFIED,
				IOobject::NO_WRITE
			)
		)
	);

	// Construct dielectric and viscocity values from dictionaries
	//(CHANGED TO NON-CONSTANT ALLOW VARYING VALUES BASED ON CONCENTRATION,ETC!!!!)
	epsilonR_Fluid.set(i, new dimensionedScalar(FluidElectroKineticProperties[i].lookup("epsilonR")));
	epsilon_Fluid.set(i,new dimensionedScalar(epsilonR_Fluid[i] * epsilonZero_));
	mu_Fluid.set(i,new dimensionedScalar(FluidElectroKineticProperties[i].lookup("mu")));

}

Info<< "Creating fluid/solid fields\n" << endl;

#include "createFluidFields.H"
#include "createSolidFields.H"
#include "writeStats.H"

const dimensionedScalar unitDensity
(
    "rho",
    dimensionSet(1,-3,0,0,0,0,0),
    scalar(1)
);

// Storage for max Courant # of all regions and species (used for DTS/adjustTimeStep)
scalar CoNum = 0.0;
