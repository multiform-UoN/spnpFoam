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

Class
    Foam::newMappedPatchFieldBase

Description
    Improvement of mappedPatchField to extend average functionality

    Example usage:

    {
        fieldName           T;          // default is same as fvPatchField
        setAverage          "scale";
        average             1.0;        // only if setAverage=true
        interpolationScheme cellPoint;  // default is cell
    }

SourceFiles
    newMappedPatchFieldBase.C

\*---------------------------------------------------------------------------*/

#ifndef newMappedPatchFieldBase_H
#define newMappedPatchFieldBase_H

#include "fixedValueFvPatchFields.H"
#include "volFieldsFwd.H"
#include "mappedPatchBase.H"
#include "interpolationCell.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mappedPatchBase;
template<class> class interpolation;

/*---------------------------------------------------------------------------*\
                  Class newMappedPatchFieldBase Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class newMappedPatchFieldBase
{

protected:

    // Protected data

        //- Mapping engine
        const mappedPatchBase& mapper_;

        //- Underlying patch field
        const fvPatchField<Type>& patchField_;

        //- Name of field to sample
        word fieldName_;

        //- "scale", "shift", "none"
        const word setAverage_;

        //- Average value the newMapped field is adjusted to maintain if
        //  setAverage_ is set true
        const Type average_;

        //- Interpolation scheme to use for nearestcell mode
        word interpolationScheme_;

public:

    // Constructors

        //- Construct from components
        newMappedPatchFieldBase
        (
            const mappedPatchBase& mapper,
            const fvPatchField<Type>& patchField,
            const word& fieldName,
            const word setAverage,
            const Type average,
            const word& interpolationScheme
        );

        //- Construct from dictionary
        newMappedPatchFieldBase
        (
            const mappedPatchBase& mapper,
            const fvPatchField<Type>& patchField,
            const dictionary& dict
        );

        //- Construct empty
        newMappedPatchFieldBase
        (
            const mappedPatchBase& mapper,
            const fvPatchField<Type>& patchField
        );

        //- Construct copy
        newMappedPatchFieldBase
        (
            const newMappedPatchFieldBase<Type>& mapper
        );

        //- Construct copy, resetting patch and field
        newMappedPatchFieldBase
        (
            const mappedPatchBase& mapper,
            const fvPatchField<Type>& patchField,
            const newMappedPatchFieldBase<Type>& base
        );


    //- Destructor
    virtual ~newMappedPatchFieldBase<Type>()
    {}


    // Member Functions

        //- Field to sample. Either on my or nbr mesh
        const GeometricField<Type, fvPatchField, volMesh>& sampleField() const;
        //- Sample generic fields
        template<class T>
        const GeometricField<T, fvPatchField, volMesh>&
        sampleVolField(word& name) const;

        template<class T>
        const GeometricField<T, fvsPatchField, surfaceMesh>&
        sampleSurField(word& name) const;


        //- Value of nbr patch
        virtual tmp<Field<Type>> newMappedField() const;

        //- Gradient on nbr patch
        virtual tmp<Field<Type>> newMappedGrad() const;

        //- detaCoeffs on nbr patch
        virtual tmp<Field<scalar>> newMappedDelta() const;

        //- internalField on nbr patch
        virtual tmp<Field<Type>> newMappedInternalField() const;

        //- get vol field on nbr patch
        template<class T>
        tmp<Field<T>> mappedVolField(word& name) const;

        //- get gradient of vol field on nbr patch
        template<class T>
        tmp<Field<T>> mappedGradField(word& name) const;

        //- get internal vol field on nbr patch
        template<class T>
        tmp<Field<T>> mappedInternalField(word& name) const;

        //- get surface field on nbr patch
        template<class T>
        tmp<Field<T>> mappedSurField(word& name) const;

        //- Write
        virtual void write(Ostream&) const;
};

#include "newMappedPatchFieldBaseTemplates.C"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "newMappedPatchFieldBase.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
