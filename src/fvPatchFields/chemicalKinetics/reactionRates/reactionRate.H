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


Class
    Foam::reactionRates::reactionRate

Description
    Base abstract class for reaction rates

Authors:
    Federico Municchi, Robert Barnett, Nottingham (2020)
\*---------------------------------------------------------------------------*/

#ifndef REACTION_RATE_H
#define REACTION_RATE_H

#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "fvCFD.H"

namespace Foam
{
namespace reactionRates
{
class reactionRate
{
private:

    //- Dictionary with settings
    dictionary  dict_;

public:

    //- Runtime type information
    TypeName("reactionRate");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        reactionRate,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );

    //- Constructor
        reactionRate( const dictionary& dict );

    //- Construct and return a clone
        autoPtr<reactionRate> clone() const
        {
            autoPtr<reactionRate> temp(New(this->dict_));
            return temp;
        }

    //- Destructor
        virtual ~reactionRate();

    //- Selector
        static autoPtr<reactionRate> New
        (
            const dictionary& dict
        );

    //- Member functions

        //- Calculate reaction rate
        virtual tmp<scalarField> value(const scalarField& S,const scalarField& F) = 0;

        //- Return derivative w.r.t. cS (solid side concentration field)
        virtual tmp<scalarField> ddcS(const scalarField& S,const scalarField& F) = 0;

        //- Return derivative w.r.t. cF (fluid side concentration field)
        virtual tmp<scalarField> ddcF(const scalarField& S,const scalarField& F) = 0;

        //- Return dictionary
        dictionary& dict() {return dict_;};
	
	//- Write 
	virtual void write(Ostream& o) const;

};
//--------------- End namespace Foam ------------------------------------------/
}
//----------------End namespace retentionRates --------------------------------/
}

#endif