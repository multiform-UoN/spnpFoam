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


\*---------------------------------------------------------------------------*/

#include "RobinPhiFvPatchScalarField.H"
// #include "fvPatchFieldMapper.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RobinPhiFvPatchScalarField::RobinPhiFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    RobinFvPatchScalarField(p, iF),
    phiName_("phi"),
    RobinKeff_(p.size())
{}


Foam::RobinPhiFvPatchScalarField::RobinPhiFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    RobinFvPatchScalarField(p, iF, dict),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    RobinKeff_(p.size())
{}


Foam::RobinPhiFvPatchScalarField::RobinPhiFvPatchScalarField
(
    const RobinPhiFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    RobinFvPatchScalarField(p, iF),
    phiName_(ptf.phiName_),
    RobinKeff_(mapper(ptf.RobinKeff_))
{}


Foam::RobinPhiFvPatchScalarField::RobinPhiFvPatchScalarField
(
    const RobinPhiFvPatchScalarField& ptf
)
:
    RobinFvPatchScalarField(ptf),
    phiName_(ptf.phiName_),
    RobinKeff_(ptf.RobinKeff_)
{}


Foam::RobinPhiFvPatchScalarField::RobinPhiFvPatchScalarField
(
    const RobinPhiFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    RobinFvPatchScalarField(ptf, iF),
    phiName_(ptf.phiName_),
    RobinKeff_(ptf.RobinKeff_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::RobinPhiFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper&   m
)
{
    RobinFvPatchScalarField::autoMap(m);
    m(RobinKeff_,RobinKeff_);
}

void Foam::RobinPhiFvPatchScalarField::rmap
(
    const fvPatchField<scalar>& ptf,
    const labelList& addr
)
{
    RobinFvPatchScalarField::rmap(ptf,addr);

    const RobinPhiFvPatchScalarField& mptf =
        refCast<const RobinPhiFvPatchScalarField>(ptf);

    RobinKeff_.rmap(mptf.RobinKeff_,addr);
}

void Foam::RobinPhiFvPatchScalarField::write(Ostream& os) const
{
    RobinFvPatchScalarField::write(os);
    writeEntry(os, "phi", phiName_);
    writeEntry(os, "RobinKeff", RobinKeff_);
}

// void Foam::RobinPhiFvPatchScalarField::evaluate
// (
//     const Pstream::commsTypes commsType
// )
// {
//
//     const scalarField& RobinK = RobinFvPatchScalarField::RobinK();
//
//     const fvsPatchField<scalar>& phip =
//         patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);
//
//     //- Calculate effective Robin coefficient
//     RobinKeff_ = phip/patch().magSf()
//                   + RobinK;
//
//     //- Evaluate Robin boundary condition
//     RobinFvPatchScalarField::evaluate();
//
// }


void Foam::RobinPhiFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
      return;
    }

    const scalarField& RobinK = RobinFvPatchScalarField::RobinK();

    const fvsPatchField<scalar>& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    //- Calculate effective Robin coefficient
    RobinKeff_ = phip/patch().magSf()
                  + RobinK;

    //- Evaluate Robin boundary condition
    RobinFvPatchScalarField::updateCoeffs();

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        RobinPhiFvPatchScalarField
    );
}

// ************************************************************************* //
