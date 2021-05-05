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
#include "linear.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace reactionRates
{
    defineTypeNameAndDebug(linear, 0);

    addToRunTimeSelectionTable
    (
        reactionRate,
        linear,
        dictionary
    );
}
}

using namespace Foam;

/*------------------------  Constructors  ------------------------------------*/
reactionRates::linear::linear
(
    const dictionary& dict
)
:
    reactionRate(dict),
    fR_(readScalar(dict.lookup("forward_rate"))),
    bR_(readScalar(dict.lookup("backward_rate")))
{}
/*------------------------  Destructors  -------------------------------------*/
reactionRates::linear::~linear()
{}
/*------------------------  Member functions ---------------------------------*/
tmp<scalarField>
reactionRates::linear::value(const scalarField& S, const scalarField& F)
{
    tmp<scalarField> tf (

        new scalarField(S.size())

    );

    scalarField& f = tf.ref();

    f = fR_*S - bR_*F;

    return tf;
}

tmp<scalarField>
reactionRates::linear::ddcF(const scalarField& S, const scalarField& F)
{
    tmp<scalarField> tf (

        new scalarField(S.size(),-bR_)

    );

    return tf;
}

tmp<scalarField>
reactionRates::linear::ddcS(const scalarField& S, const scalarField& F)
{
    tmp<scalarField> tf (

        new scalarField(S.size(),fR_)
    );

    return tf;
}