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

internalField   uniform 0;

boundaryField
{
    inlet
    {
      type            fixedValue;
      value           uniform 0;
    }
    outlet
    {
      type            fixedValue;
      value           uniform 0;
    }
    lowerWall
    {
      type            symmetry
    }
    upperWall
    {
      type            symmetry
    }
    emptyWall
    {
      type            empty;
    }
    solid_to_fluid
    {
      type            fixedValue;
      value           uniform 0;
    }

    #includeEtc "caseDicts/setConstraintTypes"
}


// ************************************************************************* //
