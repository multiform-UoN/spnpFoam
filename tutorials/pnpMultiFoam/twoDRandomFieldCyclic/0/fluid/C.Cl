/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  6
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      C.Cl;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#include "$FOAM_CASE/constant/initialConditions";

dimensions      [ 0 -3 0 0 1 0 0 ];

internalField   uniform $cFluid;

boundaryField
{
    inlet
    {
        type            cyclic;
    }
    outlet
    {
        type            cyclic;
    }
    lowerWall
    {
        type            cyclic;
    }
    upperWall
    {
        type            cyclic;
    }
    emptyWall
    {
        type            empty;
    }
    fluid_to_solid
    {
        type            noFlux;
        value           $internalField;
    }

    #includeEtc "caseDicts/setConstraintTypes"
}


// ************************************************************************* //
