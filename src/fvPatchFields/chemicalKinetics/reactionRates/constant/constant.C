/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2019
     \\/     M anipulation  | Matteo Icardi, Federico Municchi
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.


Authors:
    Federico Municchi, Robert Barnett Nottingham (2020)
\*---------------------------------------------------------------------------*/
#include "constant.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace reactionRates
{
    defineTypeNameAndDebug(constant, 0);

    addToRunTimeSelectionTable
    (
        reactionRate,
        constant,
        dictionary
    );
}
}

using namespace Foam;

/*------------------------  Constructors  ------------------------------------*/
reactionRates::constant::constant
(
    const dictionary& dict
)
:
    reactionRate(dict)
{}
/*------------------------  Destructors  -------------------------------------*/
reactionRates::constant::~constant()
{}
/*------------------------  Member functions ---------------------------------*/
tmp<scalarField>
reactionRates::constant::value(const scalarField& S, const scalarField& F)
{
    tmp<scalarField> tf (

        new scalarField(S.size(), scalar(1))

    );

    return tf;
}

tmp<scalarField>
reactionRates::constant::ddcF(const scalarField& S, const scalarField& F)
{
    tmp<scalarField> tf (

        new scalarField(S.size(),scalar(0))

    );

    return tf;
}

tmp<scalarField>
reactionRates::constant::ddcS(const scalarField& S, const scalarField& F)
{
    tmp<scalarField> tf (

        new scalarField(S.size(),scalar(0))
    );

    return tf;
}