/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2019 OpenFOAM Foundation
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

#include "newMappedPatchFieldBase.H"
#include "mappedPatchBase.H"
#include "interpolationCell.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
newMappedPatchFieldBase<Type>::newMappedPatchFieldBase
(
    const mappedPatchBase& mapper,
    const fvPatchField<Type>& patchField,
    const word& fieldName,
    const word setAverage,
    const Type average,
    const word& interpolationScheme
)
:
    mapper_(mapper),
    patchField_(patchField),
    fieldName_(fieldName),
    setAverage_(setAverage),
    average_(average),
    interpolationScheme_(interpolationScheme)
{}


template<class Type>
newMappedPatchFieldBase<Type>::newMappedPatchFieldBase
(
    const mappedPatchBase& mapper,
    const fvPatchField<Type>& patchField,
    const dictionary& dict
)
:
    mapper_(mapper),
    patchField_(patchField),
    fieldName_
    (
        dict.template lookupOrDefault<word>
        (
            "field",
            patchField_.internalField().name()
        )
    ),
    setAverage_(dict.lookupOrDefault<word>("setAverage","none")),
    average_(dict.lookupOrDefault<Type>("average",pTraits<Type>::one)),
    interpolationScheme_(interpolationCell<Type>::typeName)
{
    if (mapper_.mode() == mappedPatchBase::NEARESTCELL)
    {
        dict.lookup("interpolationScheme") >> interpolationScheme_;
    }
}


template<class Type>
newMappedPatchFieldBase<Type>::newMappedPatchFieldBase
(
    const mappedPatchBase& mapper,
    const fvPatchField<Type>& patchField
)
:
    mapper_(mapper),
    patchField_(patchField),
    fieldName_(patchField_.internalField().name()),
    setAverage_("none"),
    average_(Zero),
    interpolationScheme_(interpolationCell<Type>::typeName)
{}


template<class Type>
newMappedPatchFieldBase<Type>::newMappedPatchFieldBase
(
    const newMappedPatchFieldBase<Type>& mapper
)
:
    mapper_(mapper.mapper_),
    patchField_(mapper.patchField_),
    fieldName_(mapper.fieldName_),
    setAverage_(mapper.setAverage_),
    average_(mapper.average_),
    interpolationScheme_(mapper.interpolationScheme_)
{}


template<class Type>
newMappedPatchFieldBase<Type>::newMappedPatchFieldBase
(
    const mappedPatchBase& mapper,
    const fvPatchField<Type>& patchField,
    const newMappedPatchFieldBase<Type>& base
)
:
    mapper_(mapper),
    patchField_(patchField),
    fieldName_(base.fieldName_),
    setAverage_(base.setAverage_),
    average_(base.average_),
    interpolationScheme_(base.interpolationScheme_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const GeometricField<Type, fvPatchField, volMesh>&
newMappedPatchFieldBase<Type>::sampleField() const
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    const fvMesh& nbrMesh = refCast<const fvMesh>(mapper_.sampleMesh());

    if (mapper_.sameRegion())
    {
        if (fieldName_ == patchField_.internalField().name())
        {
            // Optimisation: bypass field lookup
            return
                dynamic_cast<const fieldType&>
                (
                    patchField_.internalField()
                );
        }
        else
        {
            const fvMesh& thisMesh = patchField_.patch().boundaryMesh().mesh();
            return thisMesh.template lookupObject<fieldType>(fieldName_);
        }
    }
    else
    {
        return nbrMesh.template lookupObject<fieldType>(fieldName_);
    }
}


template<class Type>
tmp<Field<Type>> newMappedPatchFieldBase<Type>::newMappedField() const
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag + 1;

//    const fvMesh& thisMesh = patchField_.patch().boundaryMesh().mesh();
    const fvMesh& nbrMesh = refCast<const fvMesh>(mapper_.sampleMesh());

    // Result of obtaining remote values
    tmp<Field<Type>> tnewValues(new Field<Type>(0));
    Field<Type>& newValues = tnewValues.ref();

    if (
      (mapper_.mode() == mappedPatchBase::NEARESTPATCHFACE)
      ||
      (mapper_.mode() == mappedPatchBase::NEARESTPATCHFACEAMI)
      )
        {
            const label nbrPatchID =
                nbrMesh.boundaryMesh().findPatchID(mapper_.samplePatch());

            if (nbrPatchID < 0)
            {
                FatalErrorInFunction
                 << "Unable to find sample patch " << mapper_.samplePatch()
                 << " in region " << mapper_.sampleRegion()
                 << " for patch " << patchField_.patch().name() << nl
                 << abort(FatalError);
            }

            const fieldType& nbrField = sampleField();

            newValues = nbrField.boundaryField()[nbrPatchID];
            mapper_.distribute(newValues);

        }
        else
        {
            FatalErrorInFunction
             << "Use nearestPatchFace only instead of " << mapper_.mode()
             << nl << abort(FatalError);
        }

    if (setAverage_=="scale")
    {
        Type averagePsi =
            gSum(patchField_.patch().magSf()*newValues)
           /gSum(patchField_.patch().magSf());

        newValues *= mag(average_)/mag(averagePsi);
    }
    else if (setAverage_=="shift")
    {
        Type averagePsi =
            gSum(patchField_.patch().magSf()*newValues)
           /gSum(patchField_.patch().magSf());

        newValues += (average_ - averagePsi);
    }

    // Restore tag
    UPstream::msgType() = oldTag;

    return tnewValues;
}

template<class Type>
tmp<Field<Type>> newMappedPatchFieldBase<Type>::newMappedGrad() const
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag + 1;

    //const fvMesh& thisMesh = patchField_.patch().boundaryMesh().mesh();
    const fvMesh& nbrMesh = refCast<const fvMesh>(mapper_.sampleMesh());

    // Result of obtaining remote values
    tmp<Field<Type>> tnewValues(new Field<Type>(0));
    Field<Type>& newValues = tnewValues.ref();

    if (
      (mapper_.mode() == mappedPatchBase::NEARESTPATCHFACE)
      ||
      (mapper_.mode() == mappedPatchBase::NEARESTPATCHFACEAMI)
      )
        {
            const label nbrPatchID =
                nbrMesh.boundaryMesh().findPatchID(mapper_.samplePatch());

            if (nbrPatchID < 0)
            {
                FatalErrorInFunction
                 << "Unable to find sample patch " << mapper_.samplePatch()
                 << " in region " << mapper_.sampleRegion()
                 << " for patch " << patchField_.patch().name() << nl
                 << abort(FatalError);
            }

            const fieldType& nbrField = sampleField();

            newValues =
              nbrField.boundaryField()[nbrPatchID].patch().deltaCoeffs() *
              (
              nbrField.boundaryField()[nbrPatchID] -
              nbrField.boundaryField()[nbrPatchID].patchInternalField()
              );
            mapper_.distribute(newValues);

        }
        else
        {
            FatalErrorInFunction
             << "Use nearestPatchFace only instead of " << mapper_.mode()
             << nl << abort(FatalError);
        }


    if (setAverage_=="scale")
    {
        Type averagePsi =
            gSum(patchField_.patch().magSf()*newValues)
           /gSum(patchField_.patch().magSf());

        newValues *= mag(average_)/mag(averagePsi);
    }
    else if (setAverage_=="shift")
    {
        Type averagePsi =
            gSum(patchField_.patch().magSf()*newValues)
           /gSum(patchField_.patch().magSf());

        newValues += (average_ - averagePsi);
    }

    // Restore tag
    UPstream::msgType() = oldTag;

    return tnewValues;
}


template<class Type>
tmp<Field<scalar>> newMappedPatchFieldBase<Type>::newMappedDelta() const
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag + 1;

    // Result of obtaining remote values
    tmp<Field<scalar>> tnewValues(new Field<scalar>(0));
    Field<scalar>& newValues = tnewValues.ref();

    const fvMesh& nbrMesh = refCast<const fvMesh>(mapper_.sampleMesh());

    if (
      (mapper_.mode() == mappedPatchBase::NEARESTPATCHFACE)
      ||
      (mapper_.mode() == mappedPatchBase::NEARESTPATCHFACEAMI)
      )
        {
            const label nbrPatchID =
                nbrMesh.boundaryMesh().findPatchID(mapper_.samplePatch());

            if (nbrPatchID < 0)
            {
                FatalErrorInFunction
                 << "Unable to find sample patch " << mapper_.samplePatch()
                 << " in region " << mapper_.sampleRegion()
                 << " for patch " << patchField_.patch().name() << nl
                 << abort(FatalError);
            }

            const fieldType& nbrField = sampleField();

            newValues =
              nbrField.boundaryField()[nbrPatchID].patch().deltaCoeffs();
            mapper_.distribute(newValues);

        }
        else
        {
            FatalErrorInFunction
             << "Use nearestPatchFace only instead of " << mapper_.mode()
             << nl << abort(FatalError);
        }

        // Restore tag
        UPstream::msgType() = oldTag;

        return tnewValues;
}

template<class Type>
tmp<Field<Type>> newMappedPatchFieldBase<Type>::newMappedInternalField() const
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag + 1;

    // Result of obtaining remote values
    tmp<Field<Type>> tnewValues(new Field<Type>(0));
    Field<Type>& newValues = tnewValues.ref();

    const fvMesh& nbrMesh = refCast<const fvMesh>(mapper_.sampleMesh());

    if (
      (mapper_.mode() == mappedPatchBase::NEARESTPATCHFACE)
      ||
      (mapper_.mode() == mappedPatchBase::NEARESTPATCHFACEAMI)
      )
        {
            const label nbrPatchID =
                nbrMesh.boundaryMesh().findPatchID(mapper_.samplePatch());

            if (nbrPatchID < 0)
            {
                FatalErrorInFunction
                 << "Unable to find sample patch " << mapper_.samplePatch()
                 << " in region " << mapper_.sampleRegion()
                 << " for patch " << patchField_.patch().name() << nl
                 << abort(FatalError);
            }

            const fieldType& nbrField = sampleField();

            newValues =
              nbrField.boundaryField()[nbrPatchID].patchInternalField();
            mapper_.distribute(newValues);

        }
        else
        {
            FatalErrorInFunction
             << "Use nearestPatchFace only instead of " << mapper_.mode()
             << nl << abort(FatalError);
        }

        // Restore tag
        UPstream::msgType() = oldTag;

        return tnewValues;
}

template<class Type>
void newMappedPatchFieldBase<Type>::write(Ostream& os) const
{
    writeEntry(os, "field", fieldName_);
    writeEntry(os, "setAverage", setAverage_);
    writeEntry(os, "average", average_);
    writeEntry(os, "interpolationScheme", interpolationScheme_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
