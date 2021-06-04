/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is a derivative of OpenFOAM.

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
template<class Type>
template<class T>
const GeometricField<T, fvPatchField, volMesh>&
newMappedPatchFieldBase<Type>::sampleVolField(word& name) const
{
    typedef GeometricField<T, fvPatchField, volMesh> fieldType;

    const fvMesh& nbrMesh = refCast<const fvMesh>(mapper_.sampleMesh());

    if (mapper_.sameRegion())
    {
        const fvMesh& thisMesh = patchField_.patch().boundaryMesh().mesh();
        return thisMesh.template lookupObject<fieldType>(name);
    }
    else
    {
        return nbrMesh.template lookupObject<fieldType>(name);
    }
}

template<class Type>
template<class T>
const GeometricField<T, fvsPatchField, surfaceMesh>&
newMappedPatchFieldBase<Type>::sampleSurField(word& name) const
{
    typedef GeometricField<T, fvsPatchField, surfaceMesh> fieldType;

    const fvMesh& nbrMesh = refCast<const fvMesh>(mapper_.sampleMesh());

    if (mapper_.sameRegion())
    {
        const fvMesh& thisMesh = patchField_.patch().boundaryMesh().mesh();
        return thisMesh.template lookupObject<fieldType>(name);

    }
    else
    {
        return nbrMesh.template lookupObject<fieldType>(name);
    }
}


template<class Type>
template<class T>
tmp<Field<T>>
newMappedPatchFieldBase<Type>::mappedVolField(word& name) const
{
    typedef GeometricField<T, fvPatchField, volMesh> fieldType;

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag + 1;


    const fvMesh& nbrMesh = refCast<const fvMesh>(mapper_.sampleMesh());

    // Result of obtaining remote values
    tmp<Field<T>> tnewValues(new Field<T>(0));
    Field<T>& newValues = tnewValues.ref();

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

            const fieldType& nbrField = sampleVolField<T>(name);

            newValues = nbrField.boundaryField()[nbrPatchID];
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
template<class T>
tmp<Field<T>>
newMappedPatchFieldBase<Type>::mappedGradField(word& name) const
{
    typedef GeometricField<T, fvPatchField, volMesh> fieldType;

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

            const fieldType& nbrField = sampleVolField<T>(name);

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


    // Restore tag
    UPstream::msgType() = oldTag;

    return tnewValues;
}

template<class Type>
template<class T>
tmp<Field<T>>
newMappedPatchFieldBase<Type>::mappedInternalField(word& name) const
{
    typedef GeometricField<T, fvPatchField, volMesh> fieldType;

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

            const fieldType& nbrField = sampleVolField<T>(name);

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
template<class T>
tmp<Field<T>>
newMappedPatchFieldBase<Type>::mappedSurField(word& name) const
{
    typedef GeometricField<T, fvsPatchField, surfaceMesh> fieldType;

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag + 1;

    const fvMesh& nbrMesh = refCast<const fvMesh>(mapper_.sampleMesh());


    // Result of obtaining remote values
    tmp<Field<T>> tnewValues(new Field<T>(0));
    Field<T>& newValues = tnewValues.ref();

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

            const fieldType& nbrField = sampleSurField<T>(name);

            newValues = nbrField.boundaryField()[nbrPatchID];
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
