//Run through all species and update Courant numbers
forAll(species,sp){

	//Gather pointers to specific ion species sp
                                                                    //volScalarField&   C(cPtrL[sp]);
    dimensionedScalar D(species[sp].lookup("D"));
    dimensionedScalar Z(species[sp].lookup("Z"));

    scalarField sumPhi
    (
        fvc::surfaceSum(mag(phi - fvc::flux(e_*fvc::grad(V)/(k_*T))*Z*D))().primitiveField()
    );

    // Species Maximum Co #
    speciesCoNum[sp] = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();

    // Species Mean Co #
    speciesMeanCoNum[sp] = 0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();

    //Set CoNum to the maximum Co # of all the species
    if(CoNum<speciesCoNum[sp])
    {
    	CoNum = speciesCoNum[sp];
    }

    Info<< "Species "<< specEntries[sp].keyword() << " Courant Number mean: " << speciesMeanCoNum[sp]
        << " max: " << speciesCoNum[sp] << endl;
}
