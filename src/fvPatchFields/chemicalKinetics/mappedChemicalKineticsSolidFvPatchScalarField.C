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

#include "mappedChemicalKineticsSolidFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"


#define MAX_NEWTON_ITER 1000
#define NEWTON_TOLERANCE 1e-14
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mappedChemicalKineticsSolidFvPatchScalarField::mappedChemicalKineticsSolidFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    RobinFvPatchScalarField(p, iF),
    newMappedPatchFieldBase<scalar>(this->mapper(p, iF), *this),
    localPhi_(""),
    mappedPhi_(""),
    RobinKeff_(p.size()),
    RobinFeff_(p.size()),
    Dfluid_(p.size()),
    R_(NULL)
{}


Foam::mappedChemicalKineticsSolidFvPatchScalarField::mappedChemicalKineticsSolidFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    RobinFvPatchScalarField(p, iF, dict),
    newMappedPatchFieldBase<scalar>(this->mapper(p, iF), *this, dict),
    localPhi_(dict.lookupOrDefault<word>("localPhi", "")),
    mappedPhi_(dict.lookupOrDefault<word>("mappedPhi","")), // CHANGE THESE SO THAT THEY DEFAULT TO ZERO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    RobinKeff_(p.size()),
    RobinFeff_(p.size()),
    Dfluid_(p.size()),
    R_
    (
        reactionRates::reactionRate::New
        (
            dict.subDict("reaction")
        )
    )
{
  if (dict.found("Dfluid"))
  {
    Dfluid_ = scalarField("Dfluid",dict,p.size());
  }
  else
  {
    Dfluid_ = scalarField(p.size(),1.0);
  }

  if (dict.found("Dsolid"))
  {
    RobinD(scalarField("Dsolid",dict,p.size()));
  }
}


Foam::mappedChemicalKineticsSolidFvPatchScalarField::mappedChemicalKineticsSolidFvPatchScalarField
(
    const mappedChemicalKineticsSolidFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    RobinFvPatchScalarField(p, iF),
    newMappedPatchFieldBase<scalar>(this->mapper(p, iF), *this, ptf),
    localPhi_(ptf.localPhi_),
    mappedPhi_(ptf.mappedPhi_),
    RobinKeff_(mapper(ptf.RobinKeff_)),
    RobinFeff_(mapper(ptf.RobinFeff_)),
    Dfluid_(mapper(ptf.Dfluid_)),
    R_(ptf.R_,false)
{}


Foam::mappedChemicalKineticsSolidFvPatchScalarField::mappedChemicalKineticsSolidFvPatchScalarField
(
    const mappedChemicalKineticsSolidFvPatchScalarField& ptf
)
:
    RobinFvPatchScalarField(ptf),
    newMappedPatchFieldBase<scalar>(ptf),
    localPhi_(ptf.localPhi_),
    mappedPhi_(ptf.mappedPhi_),
    RobinKeff_(ptf.RobinKeff_),
    RobinFeff_(ptf.RobinFeff_),
    Dfluid_(ptf.Dfluid_),
    R_(ptf.R_,false)
{}


Foam::mappedChemicalKineticsSolidFvPatchScalarField::mappedChemicalKineticsSolidFvPatchScalarField
(
    const mappedChemicalKineticsSolidFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    RobinFvPatchScalarField(ptf, iF),
    newMappedPatchFieldBase<scalar>(this->mapper(this->patch(), iF), *this, ptf),
    localPhi_(ptf.localPhi_),
    mappedPhi_(ptf.mappedPhi_),
    RobinKeff_(ptf.RobinKeff_),
    RobinFeff_(ptf.RobinFeff_),
    Dfluid_(ptf.Dfluid_),
    R_(ptf.R_,false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::mappedPatchBase& Foam::mappedChemicalKineticsSolidFvPatchScalarField::mapper
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

void Foam::mappedChemicalKineticsSolidFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper&   m
)
{
    RobinFvPatchScalarField::autoMap(m);
    //m(localPhi_,localPhi_);
    //m(mappedPhi_,mappedPhi_),
    m(RobinKeff_,RobinKeff_);
    m(RobinFeff_,RobinFeff_);
    //m(R_,R_);
    m(Dfluid_,Dfluid_);
}

void Foam::mappedChemicalKineticsSolidFvPatchScalarField::rmap
(
    const fvPatchField<scalar>& ptf,
    const labelList& addr
)
{
    RobinFvPatchScalarField::rmap(ptf,addr);

    const mappedChemicalKineticsSolidFvPatchScalarField& mptf =
        refCast<const mappedChemicalKineticsSolidFvPatchScalarField>(ptf);

    //localPhi_.rmap(mptf.localPhi_,addr);
    //mappedPhi_.rmap(mptf.mappedPhi_,addr);
    RobinKeff_.rmap(mptf.RobinKeff_,addr);
    RobinFeff_.rmap(mptf.RobinFeff_,addr);
    R_ = mptf.R_;
    Dfluid_.rmap(mptf.Dfluid_,addr);
}

void Foam::mappedChemicalKineticsSolidFvPatchScalarField::write(Ostream& os) const
{
    RobinFvPatchScalarField::write(os);
    newMappedPatchFieldBase::write(os);
    writeEntry(os, "localPhi_", localPhi_);
    writeEntry(os,"mappedPhi_",mappedPhi_);
//    writeEntry(os, "RobinKeff", RobinKeff_);
//    writeEntry(os, "RobinFeff", RobinFeff_);
    writeEntry(os, "Dfluid", Dfluid_);
    // ADD ALSO THE REACTION
}


void Foam::mappedChemicalKineticsSolidFvPatchScalarField::updateCoeffs()
{
//    Info << "update " << this->updated() <<  endl;

    if (this->updated())
    {
        return;
    }

    // Set the update flag to true
    //RobinFvPatchScalarField::updateCoeffs();

    // Reference to internal mesh of this patch
    const fvMesh& thisMesh = this->internalField().mesh();
    // Reference to internal mesh of mapped (other) patch
    const fvMesh& nbrMesh = refCast<const fvMesh>(mapper_.sampleMesh());

    const scalarField& RobinK0 = RobinFvPatchScalarField::RobinK();
    const scalarField& RobinF0 = RobinFvPatchScalarField::RobinF();

    // Storage of local phi divided by face magnitudes (stays zero if localPhi_ not defined/non-existant)
    scalarField lPhiNorm(patch().size(),scalar(0));

    // Storage of mapped phi divided by face magnitudes (stays zero if mappedPhi_ not defined/non-existant)
    scalarField mPhiNorm(patch().size(),scalar(0));

    //Look up flux variable on this patch if localPhi_ non-empty and found
    if(localPhi_ != "" && thisMesh.objectRegistry::template foundObject<surfaceScalarField>(localPhi_)){

        lPhiNorm = patch().lookupPatchField<surfaceScalarField, scalar>(localPhi_)/patch().magSf();
    }

    // Look up flux variable on mapped patch if mappedPhi_ non-empty and found
    if(mappedPhi_ != "" && nbrMesh.objectRegistry::template foundObject<surfaceScalarField>(mappedPhi_)){

        mPhiNorm = mappedSurField<scalar>(mappedPhi_)/patch().magSf();
    }

    // Store field value on nbr patch
    tmp<Field<scalar>> CFluid = this->newMappedField();
    // Create reference to field value in tmp wrapper
    scalarField& CFluidref(CFluid.ref());

    // Iterator for Newton method
    label n = 0;

    // Current and previous values of field
    scalarField& C(*this);
    scalarField C0(*this);

    do
    {
        if (n>MAX_NEWTON_ITER)
        {
            FatalErrorInFunction
            << "Max number of Newton iterations reached for "
            << "mappedChemicalKineticsFvPatchScalarField: "
            << this->patch().name()
            << exit(FatalError);
        }

        // Store field value of previous iteration
        C0 = C;

        // Compute next guess of effective Robin coefficients
        RobinKeff_ =    (
                            R_().ddcS(C,CFluidref)*R_().ddcF(C,CFluidref)
                        )
                        /
                        (
                            mPhiNorm - Dfluid_*this->newMappedDelta() + R_().ddcF(C,CFluidref)
                        )

                        - R_().ddcS(C,CFluidref) + lPhiNorm + RobinK0;


        RobinFeff_ =    - R_().value(C,CFluidref) + CFluidref*R_().ddcF(C,CFluidref) + C*R_().ddcS(C,CFluidref)

                        - R_().ddcF(C,CFluidref)*
                        (
                            - Dfluid_*this->newMappedDelta()*this->newMappedInternalField() - R_().value(C,CFluidref)

                            + CFluidref*R_().ddcF(C,CFluidref) + C*R_().ddcS(C,CFluidref)
                        )
                        /
                        (
                            mPhiNorm + R_().ddcF(C,CFluidref) - Dfluid_*this->newMappedDelta()
                        )

                        + RobinF0;

        // - Test to enforce continuity of fluxes and preserve mass
        // RobinKeff_ = RobinK0;//mPhiNorm - this->newMappedDelta()*Dfluid_ + RobinK0;
        // RobinFeff_ = mPhiNorm -Dfluid_*this->newMappedDelta()*(CFluidref - this->newMappedInternalField())
        //           + RobinF0;

        // Set updated_ flag of coeffs to true
        RobinFvPatchScalarField::updateCoeffs();

        // Evaluate patch to update to next field value iteration
        RobinFvPatchScalarField::evaluate();

        // Debug
        // Info<< "solid flux ==" << lPhiNorm -this->RobinD()*this->snGrad()<<endl;
        // Info<< "solid advective flux comp ==" << lPhiNorm << endl;
        // Info<< "solid diffusive flux comp ==" << -this->RobinD()*this->snGrad()<< endl;

        // Info<< "error == " << Foam::sqrt(gSum(pow(C0-C,2))) << endl;

        // Increment iterator
        n++;
    }while(Foam::sqrt(gSum(pow(C0-C,2)))>NEWTON_TOLERANCE);

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        mappedChemicalKineticsSolidFvPatchScalarField
    );
}

// ************************************************************************* //
