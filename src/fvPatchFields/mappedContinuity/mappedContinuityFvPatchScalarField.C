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

#include "mappedContinuityFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mappedContinuityFvPatchScalarField::mappedContinuityFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    RobinFvPatchScalarField(p, iF),
    newMappedPatchFieldBase<scalar>(this->mapper(p, iF), *this),
    RobinKeff_(p.size()),
    RobinFeff_(p.size()),
    Dnbr_(p.size())
{}


Foam::mappedContinuityFvPatchScalarField::mappedContinuityFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    RobinFvPatchScalarField(p, iF, dict),
    newMappedPatchFieldBase<scalar>(this->mapper(p, iF), *this, dict),
    RobinKeff_(p.size()),
    RobinFeff_(p.size()),
    Dnbr_("Dnbr", dict, p.size())
{
    //updateCoeffs();
}


Foam::mappedContinuityFvPatchScalarField::mappedContinuityFvPatchScalarField
(
    const mappedContinuityFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    RobinFvPatchScalarField(p, iF),
    newMappedPatchFieldBase<scalar>(this->mapper(p, iF), *this, ptf),
    RobinKeff_(mapper(ptf.RobinKeff_)),
    RobinFeff_(mapper(ptf.RobinFeff_)),
    Dnbr_(mapper(ptf.Dnbr_))
{}


Foam::mappedContinuityFvPatchScalarField::mappedContinuityFvPatchScalarField
(
    const mappedContinuityFvPatchScalarField& ptf
)
:
    RobinFvPatchScalarField(ptf),
    newMappedPatchFieldBase<scalar>(ptf),
    RobinKeff_(ptf.RobinKeff_),
    RobinFeff_(ptf.RobinFeff_),
    Dnbr_(ptf.Dnbr_)
{}


Foam::mappedContinuityFvPatchScalarField::mappedContinuityFvPatchScalarField
(
    const mappedContinuityFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    RobinFvPatchScalarField(ptf, iF),
    newMappedPatchFieldBase<scalar>(this->mapper(this->patch(), iF), *this, ptf),
    RobinKeff_(ptf.RobinKeff_),
    RobinFeff_(ptf.RobinFeff_),
    Dnbr_(ptf.Dnbr_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::mappedPatchBase& Foam::mappedContinuityFvPatchScalarField::mapper
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
{
    if (!isA<mappedPatchBase>(p.patch()))
    {
        FatalErrorInFunction
            << "' not type '" << mappedPatchBase::typeName << "'"
            << "\n    for patch " << p.patch().name()
            << " of field " << iF.name()
            << " in file " << iF.objectPath()
            << exit(FatalError);
    }
    return refCast<const mappedPatchBase>(p.patch());
}

// void Foam::mappedContinuityFvPatchScalarField::autoMap
// (
//     const fvPatchFieldMapper&   m
// )
// {
//     RobinFvPatchScalarField::autoMap(m);
// //    m(phiName_,phiName_);
// //    m(RobinKeff_,RobinKeff_);
// }
//
// void Foam::mappedContinuityFvPatchScalarField::rmap
// (
//     const fvPatchField<scalar>& ptf,
//     const labelList& addr
// )
// {
//     RobinFvPatchScalarField::rmap(ptf,addr);
//
//     const mappedContinuityFvPatchScalarField& mptf =
//         refCast<const mappedContinuityFvPatchScalarField>(ptf);
//
// //    phiName_.rmap(mptf.phiName_,addr);
// //    RobinKeff_.rmap(mptf.RobinKeff_,addr);
// }

void Foam::mappedContinuityFvPatchScalarField::write(Ostream& os) const
{
    RobinFvPatchScalarField::write(os);
    newMappedPatchFieldBase::write(os);
//    writeEntry(os, "RobinKeff", RobinKeff_);
//    writeEntry(os, "RobinFeff", RobinFeff_);
    writeEntry(os, "Dnbr", Dnbr_);
}

// void Foam::mappedContinuityFvPatchScalarField::evaluate
// (
//     const Pstream::commsTypes commsType
// )
// {
//
//     const scalarField& RobinK = RobinFvPatchScalarField::RobinK();
//     const scalarField& RobinF = RobinFvPatchScalarField::RobinF();
//
//     //- Calculate effective Robin coefficient
//     RobinKeff_ =
//                   - this->newMappedDelta()*Dnbr_
//                   + RobinK;
//     RobinFeff_ = this->newMappedDelta()*this->newMappedInternalField()*Dnbr_
//                   + RobinF;
//
//     //- Evaluate Robin boundary condition
//     RobinFvPatchScalarField::evaluate();
//
// }

void Foam::mappedContinuityFvPatchScalarField::updateCoeffs()
{
//    Info << "update " << this->updated() <<  endl;

    if (this->updated())
    {
        return;
    }

    const scalarField& RobinK = RobinFvPatchScalarField::RobinK();
    const scalarField& RobinF = RobinFvPatchScalarField::RobinF();

    //- Calculate effective Robin coefficient
    RobinKeff_ =  - this->newMappedDelta()*Dnbr_
                  + RobinK;
    RobinFeff_ = this->newMappedDelta()*this->newMappedInternalField()*Dnbr_
                  + RobinF;

    //- Evaluate Robin boundary condition
    RobinFvPatchScalarField::updateCoeffs();

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        mappedContinuityFvPatchScalarField
    );
}

// ************************************************************************* //
