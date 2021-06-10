//PtrList<volScalarField> fieldsfluid(fluidRegions.size());

PtrList<volScalarField> VFluid(fluidRegions.size());
PtrList<volScalarField> TFluid(fluidRegions.size());
PtrList<surfaceScalarField> phiNPFluid(fluidRegions.size());
PtrList<volVectorField> UFluid(fluidRegions.size());
PtrList<volScalarField> rhoFluid(fluidRegions.size());
PtrList<surfaceScalarField> phiFluid(fluidRegions.size());
PtrList<volScalarField> pFluid(fluidRegions.size());

// fvOptions per fluid region
PtrList<fv::options> fvOptionsFluid(fluidRegions.size());

// List of electro-neutrality and diffusion dominant flags
List<bool> electroNeutralityFluid(fluidRegions.size(), false);

// Lists of reference cells and values
List<label> VRefCellFluid(fluidRegions.size());
List<scalar> VRefValueFluid(fluidRegions.size());
List<label> pRefCellFluid(fluidRegions.size());
List<scalar> pRefValueFluid(fluidRegions.size());

// Pointer list towards each fluid regions concentration pointer list
PtrList<PtrList<volScalarField>> cPtrLFluid(fluidRegions.size());
PtrList<PtrList<surfaceScalarField>> phiPtrLFluid(fluidRegions.size());

// Pointer list of species dictionary list for fluid regions
PtrList<PtrList<dictionary>> speciesFluid(fluidRegions.size());

// Pointer list of species Max Co # for fluid regions
PtrList<List<scalar>> speciesCoNumFluid(fluidRegions.size());

// Pointer list of mean species Co # for fluid regions
PtrList<List<scalar>> speciesMeanCoNumFluid(fluidRegions.size());

// Pointer list of total mass
PtrList<List<scalar>> speciesTotConcFluid(fluidRegions.size());
// Volumes
List<scalar> volFluid(fluidRegions.size(),0.0);

// List a Max Co # of all species for fluid regions (initialize all to 0)
List<scalar> MaxCoFluid(fluidRegions.size(),0.0);

// Populate fluid field pointer lists
forAll(fluidRegions, i){

    // Compute volumes
    volFluid[i] = Foam::gSum(fluidRegions[i].V());
    Info<< "Fluid region " << fluidRegions[i].name()
          << " Volume " << volFluid[i] << endl;

    //Look up species from dictionary of fluid region i
    PtrList<entry> specEntries(

            FluidElectroKineticProperties[i].lookup("species")
        );

    //Set List entry to electro-neutrality flag in dictionary
    electroNeutralityFluid[i] = readBool(FluidElectroKineticProperties[i].lookup("electroNeutrality"));

    //Pre-allocate PtrList of PtrLists for species dictionaries and concentrations
    cPtrLFluid.set(
        i,
        new PtrList<volScalarField>(specEntries.size())
    );

    phiPtrLFluid.set(
        i,
        new PtrList<surfaceScalarField>(specEntries.size())
    );

    speciesFluid.set(
        i,
        new PtrList<dictionary>(specEntries.size())
    );

    Info<< "*** Reading fluid mesh properties for region "
        << fluidRegions[i].name() << nl << endl;


    //VFluid := fluid electric potential field
    Info<< "    Adding to VFluid\n" << endl;
    VFluid.set(
        i,
        new volScalarField
        (
            IOobject
            (
                "V",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ), fluidRegions[i]
        )
    );

    //pFluid := fluid pressure field
    Info<< "    Adding to pFluid\n" << endl;
    pFluid.set(
        i,
        new volScalarField
        (
            IOobject
            (
                "p",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ), fluidRegions[i]
        )
    );

    VRefCellFluid[i] = 0;
    VRefValueFluid[i] = 0.0;

    pRefCellFluid[i] = 0;
    pRefValueFluid[i] = 0.0;

    //Update RefCellFluid and RefValueFluid if need to
    setRefCell(
        VFluid[i],
        fluidRegions[i].solutionDict().subDict("PIMPLE"),
        VRefCellFluid[i],
        VRefValueFluid[i],
        fluidRegions[i].solutionDict().subDict("PIMPLE").lookupOrDefault<bool>("forceVRef",false)
    );

    setRefCell(
        pFluid[i],
        fluidRegions[i].solutionDict().subDict("PIMPLE"),
        pRefCellFluid[i],
        pRefValueFluid[i]
    );

    //TFluid := fluid temperature field
    Info<<"     Adding to TFluid\n"<< endl;
    TFluid.set(
        i,
        new volScalarField
        (
            IOobject
            (
                "T",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ), fluidRegions[i]
        )
    );



    //UFluid := fluid velocity field
    Info<<"     Adding to UFluid\n"<< endl;
    UFluid.set(
        i,
        new volVectorField
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),  fluidRegions[i]
        )
    );

    //rhoFluid := fluid charge density field
    Info<<"     Adding to rhoFluid\n"<< endl;
    rhoFluid.set(
        i,
        new volScalarField
        (
            IOobject
            (
                "rho",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i],
            dimensionedScalar(
                "rho",
                dimensionSet(0,-3,1,0,0,1,0),
                scalar(0)
            )
        )
    );

    // phiFluid := face flux field
    Info<<"     Adding to phiFluid\n"<< endl;
    phiFluid.set(
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            linearInterpolate(UFluid[i]) & fluidRegions[i].Sf()
        )
    );

    Info<<"     Adding to species list\n"<< endl;
    //- Set all species objects in list
    for(int sp=0; sp<specEntries.size(); sp++)
    {

        cPtrLFluid[i].set(
            sp,
            new volScalarField
            (
                IOobject
                (
                    specEntries[sp].keyword(),
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i]
            )

        );

        phiPtrLFluid[i].set(
            sp,
            new surfaceScalarField
            (
                IOobject
                (
                    "phi"+specEntries[sp].keyword(),
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                phiFluid[i]
            )
        );

        speciesFluid[i].set(sp,new dictionary(specEntries[sp].dict()));

    }

    // Initialize species list of Max and Mean Co #'s to zero
    speciesCoNumFluid.set(
        i,
        new List<scalar>(specEntries.size(),0.0)
        );

    speciesMeanCoNumFluid.set(
        i,
        new List<scalar>(specEntries.size(),0.0)
        );

    speciesTotConcFluid.set(
        i,
        new List<scalar>(specEntries.size(),0.0)
        );

    // Set fvOptions per fluid region (!!!!!!!!! UNCLEAR IF THIS IS CORRECT AS fv::options HAS A New FUNCTION !!!!!!!)
    Info<<"     Adding fvOptions\n"<< endl;
    fvOptionsFluid.set(
        i,
        new fv::options(fluidRegions[i])
    );

    if (!fvOptionsFluid[i].optionList::size())
    {
        Info << "No finite volume options present" << endl;
    }

    // Set the flux required for pressure per fluid region (Post-processing purposes)
    fluidRegions[i].setFluxRequired(pFluid[i].name());

}
