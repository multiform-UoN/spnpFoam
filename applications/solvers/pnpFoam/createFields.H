//- p:= Pressure field
Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

//- V := Electric potential
Info<< "Reading field V\n" << endl;
volScalarField V
(
    IOobject
    (
        "V",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//- T:= Temperature
Info<< "Reading field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


// Initialize reference value and cell
label VRefCell = 0;
scalar VRefValue = 0.0;

//Update RefCell and RefValue if need to
setRefCell
(
    V,
    mesh.solutionDict().subDict("PIMPLE"),
    VRefCell,
    VRefValue
);


//- U:= Velocity field
Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//Creates and initialises the relative face-flux field phi.
#include "createPhi.H"

Info<< "Reading/calculating Nernst-Planck face flux field phiNP\n" << endl;

surfaceScalarField phiNP
(
    IOobject
    (
        "phiNP",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    phi*scalar(0)
);


//- rho:= charge density field
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(
        "rho",
        dimensionSet(0,-3,1,0,0,1,0),
        scalar(0)
    )
);


Info<< "Reading electrokineticProperties\n" << endl;

IOdictionary electrokineticProperties
(
    IOobject
    (
        "electrokineticProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

const bool electroNeutrality(
    electrokineticProperties.lookupOrDefault<bool>("electroNeutrality",false) //readBool
);

Info<< "Electro-neutrality flag set to "<< electroNeutrality<<endl;

//- Read dynamic viscosity
const dimensionedScalar mu
(
    electrokineticProperties.lookup("mu")
);

//- Read body force parameter
const scalar bodyForce
(
    electrokineticProperties.lookupOrDefault<scalar>("bodyForce",scalar(1))
);

//- Read fluid relative dielectric
const dimensionedScalar epsilonR_
(
    electrokineticProperties.lookup("epsilonR")
);

//- Calculate effective dielectric
const dimensionedScalar epsilon_( epsilonR_ * epsilonZero_);

//- Read chemical species
PtrList<dictionary> species;

PtrList<entry> specEntries
(
    electrokineticProperties.lookup("species")
);

species.setSize(specEntries.size());

//- Create concentration fields
PtrList<volScalarField> cPtrL(species.size());

forAll(species,sp)
{

    cPtrL.set(
        sp,
        new volScalarField
        (
            IOobject
            (
                specEntries[sp].keyword(),
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )

    );

    species.set(sp,new dictionary(specEntries[sp].dict()));

}

const dimensionedScalar unitDensity
(
    "rho",
    dimensionSet(1,-3,0,0,0,0,0),
    scalar(1)
);

//Initialize list of Courant numbers for all species to 0
List<scalar> speciesCoNum(specEntries.size(),0.0);
List<scalar> speciesMeanCoNum(specEntries.size(),0.0);

//Initialize max Co # of all species to 0
scalar CoNum = 0.0;

// Set reference pressure value and location
label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("PIMPLE"), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());

#include "createMRF.H"
#include "createFvOptions.H"

// Test: functionModel::linear
//#include "testFunctionModel.H"
