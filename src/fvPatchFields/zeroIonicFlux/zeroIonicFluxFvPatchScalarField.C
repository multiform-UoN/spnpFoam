/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "zeroIonicFluxFvPatchScalarField.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvCFD.H"
#include "electrokineticConstants.H"
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zeroIonicFluxFvPatchScalarField::zeroIonicFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    zib_()
{}


Foam::zeroIonicFluxFvPatchScalarField::zeroIonicFluxFvPatchScalarField
(
    const zeroIonicFluxFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    zib_(ptf.zib_)
{}


Foam::zeroIonicFluxFvPatchScalarField::zeroIonicFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    zib_(0)
{

    if (dict.found("value") && dict.found("gradient"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
        gradient() = scalarField("gradient", dict, p.size());
    }
    else
    {
        fvPatchField<scalar>::operator=(patchInternalField());
        gradient() = 0.0;
    }

    const dictionary& elecDict =
 		db().lookupObject<IOdictionary>("electrokineticProperties");

    PtrList<entry> specEntries_(elecDict.lookup("species"));

    forAll (specEntries_, specI)
    {
       if ( specEntries_[specI].keyword() == this->internalField().name() )
        {
          dimensionedScalar zid_(specEntries_[specI].dict().lookup("Z"));
          zib_ = zid_.value();
          break;
        }
    }
}


Foam::zeroIonicFluxFvPatchScalarField::zeroIonicFluxFvPatchScalarField
(
    const zeroIonicFluxFvPatchScalarField& wbppsf
)
:
    fixedGradientFvPatchScalarField(wbppsf),
    zib_(wbppsf.zib_)
{}


Foam::zeroIonicFluxFvPatchScalarField::zeroIonicFluxFvPatchScalarField
(
    const zeroIonicFluxFvPatchScalarField& wbppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(wbppsf, iF),
    zib_(wbppsf.zib_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::zeroIonicFluxFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // const dictionary& elecDict = db().lookupObject<IOdictionary>("electrokineticProperties");

    const volScalarField&  T_(db().lookupObject<volScalarField>("T"));

    const volScalarField& V_ =
        db().lookupObject<volScalarField>("V");

	scalarField Tpatch( T_.boundaryField()[patch().index()]);
    scalarField Epatch( V_.boundaryField()[patch().index()].snGrad() );

    scalar eK_(electrokineticConstants::e_.value());
    scalar kbK_(electrokineticConstants::k_.value());


    gradient() = -(*this) * eK_ * zib_ * Epatch / (kbK_ * Tpatch);


    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::zeroIonicFluxFvPatchScalarField::evaluate(const Pstream::commsTypes)
{

	// const dictionary& elecDict = db().lookupObject<IOdictionary>("electrokineticProperties");

    const volScalarField& T_(db().lookupObject<volScalarField>("T"));

    const volScalarField& V_ =
        db().lookupObject<volScalarField>("V");

	scalarField Tpatch( T_.boundaryField()[patch().index()]);
    scalarField Epatch( V_.boundaryField()[patch().index()].snGrad() );

    scalar eK_(electrokineticConstants::e_.value());
    scalar kbK_(electrokineticConstants::k_.value());

    const fvPatchField<scalar>& psib_ = patch().lookupPatchField<volScalarField, scalar>("V");

    scalarField deltaPsi(psib_ - psib_.patchInternalField());



    Field<Foam::scalar>::operator=
    (
        this->patchInternalField() * Foam::exp( -eK_ * zib_ * deltaPsi/ (kbK_ * Tpatch) )
    );

    // This will set update_ = false and force updateCoeffs(). Gradient
    // should be updated always we call operator= in order to ensure
    // that both use the same deltaPsi.
    fvPatchField<Foam::scalar>::evaluate();

    if (!this->updated())
    {
        this->updateCoeffs();
    }

    fvPatchField<Foam::scalar>::evaluate();
}

// Render the BC fully-explicit

Foam::tmp<Field<Foam::scalar> > Foam::zeroIonicFluxFvPatchScalarField::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
     return (*this);
}

// Render the BC fully-explicit

Foam::tmp<Field<Foam::scalar> > Foam::zeroIonicFluxFvPatchScalarField::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<Field<Foam::scalar> >(new Field<Foam::scalar>(this->size(), pTraits<Foam::scalar>::zero));
}

Foam::tmp<Field<Foam::scalar> >
Foam::zeroIonicFluxFvPatchScalarField::gradientInternalCoeffs() const
{
    return tmp<Field<Foam::scalar>>
    (
        new Field<Foam::scalar>(this->size(), Zero)
    );
}

Foam::tmp<Field<Foam::scalar> >
Foam::zeroIonicFluxFvPatchScalarField::gradientBoundaryCoeffs() const
{
     return gradient();
}


void Foam::zeroIonicFluxFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        zeroIonicFluxFvPatchScalarField
    );
}

// ************************************************************************* //
