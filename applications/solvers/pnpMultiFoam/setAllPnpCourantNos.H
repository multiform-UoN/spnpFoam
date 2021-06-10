// SETS/UPDATES ALL COURANT NUMBERS OF ALL REGIONS

// Set Co #'s for all fluid regions
forAll(fluidRegions,i){

	// Set fields, mesh and options to fluid mesh i
	#include "setRegionFluidFields.H"

	forAll(species,sp)
	{

		// Gather pointers to specific ion species sp
		volScalarField&   C(cPtrL[sp]);
	    surfaceScalarField& phiNP(phiPtrL[sp]);
		dimensionedScalar D(species[sp].lookup("D"));
		dimensionedScalar Z(species[sp].lookup("Z"));

		// Update Nernst-Planck flux to fluid region i with species sp
		phiNP = phi - fvc::flux(e_*fvc::grad(V)/(k_*T))*Z*D;

		if (runTime.outputTime())
		{
			const volVectorField U("U"+species[sp].name(),fvc::reconstruct(phiNP));
			U.write();
		}

		// Take surface sum of NP flux
		scalarField sumPhi
		(
			fvc::surfaceSum(mag(phiNP))().primitiveField()
		);

		// Species Maximum Co #
		spCoNum[sp] = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();

		// Species Mean Co #
		spMeanCoNum[sp] = 0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();

		//Update MaxCoFluid to the maximum Co # of all the species in fluid region i
		if(MaxCoFluid[i]<spCoNum[sp]){

			MaxCoFluid[i] = spCoNum[sp];

		}

		Info<< "Species "<< species[sp].name()
				<<" Courant Number mean: " << spMeanCoNum[sp]
				<< " max: " << spCoNum[sp] << endl;

	}
}

// ------------ TODO
// ADD HERE THE CALCULATION OF THE DIFFUSION NUMBER


// // Set Co #'s for all solid regions
// forAll(solidRegions,i){
//
// 	// Set fields, mesh and options to solid mesh i
// 	#include "setRegionSolidFields.H"
//
// 	forAll(species,sp)
// 	{
// 		// Gather pointers to specific ion species sp
// 		volScalarField&   C(cPtrL[sp]);
// 		dimensionedScalar D(species[sp].lookup("D"));
// 		dimensionedScalar Z(species[sp].lookup("Z"));
//
// 		// Update Nernst-Planck flux to fluid region i with species sp
// 		phiNP = - fvc::flux(e_*fvc::grad(V)/(k_*T))*Z*D;
//
// 		if (runTime.outputTime())
// 		{
// 			const volVectorField U("U"+species[sp].name(),fvc::reconstruct(phiNP));
// 			U.write();
// 		}
//
// 		// Take surface sum of NP flux
// 		scalarField sumPhi
// 		(
// 			fvc::surfaceSum(mag(phiNP))().primitiveField()
// 		);
//
// 		// Species Maximum Co #
// 		spCoNum[sp] = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();
//
// 		// Species Mean Co #
// 		spMeanCoNum[sp] = 0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();
//
// 		//Update MaxCoSolid to the maximum Co # of all the species in solid region i
// 		if(MaxCoSolid[i]<spCoNum[sp])
// 		{
// 			MaxCoSolid[i] = spCoNum[sp];
// 		}
//
// 		Info<< "Species "<< species[sp].name()
// 				<<" Courant Number mean: " << spMeanCoNum[sp]
// 				<< " max: " << spCoNum[sp] << endl;
//
// 	}
// }
