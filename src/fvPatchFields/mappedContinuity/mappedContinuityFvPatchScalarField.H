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
    Foam::mappedContinuityFvPatchScalarField

Description
    This boundary condition impose continuity and continuous fluxes
    across a mapped patch for conjugate problems from the fluid side
    advective and diffusive are considered while on the nbr patch (solid)
    only diffusion is used with diffusivity Dnbr_

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    \endverbatim


SourceFiles
    mappedContinuityFvPatchScalarField.C

Authors:
    Matteo Icardi, Federico Municchi, Nottingham (2019)
\*---------------------------------------------------------------------------*/

#ifndef mappedContinuityFvPatchScalarField_H
#define mappedContinuityFvPatchScalarField_H

#include "RobinFvPatchFields.H"
#include "newMappedPatchFieldBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class mappedContinuityFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class mappedContinuityFvPatchScalarField
:
    public RobinFvPatchScalarField,
    public newMappedPatchFieldBase<scalar>
{

protected:

    //- Effective RobinK_
    scalarField RobinKeff_;

    //- Effective RobinF_
    scalarField RobinFeff_;

    //- Effective Dnbr_
    scalarField Dnbr_;

    // Protected Member Functions

        const mappedPatchBase& mapper
        (
            const fvPatch& p,
            const DimensionedField<scalar, volMesh>& iF
        );

public:

    //- Runtime type information
    TypeName("mappedContinuity");


    // Constructors

        //- Construct from patch and internal field
        mappedContinuityFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        mappedContinuityFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given mappedContinuityFvPatchScalarField
        //  onto a new patch
        mappedContinuityFvPatchScalarField
        (
            const mappedContinuityFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        mappedContinuityFvPatchScalarField
        (
            const mappedContinuityFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new mappedContinuityFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        mappedContinuityFvPatchScalarField
        (
            const mappedContinuityFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new mappedContinuityFvPatchScalarField(*this, iF)
            );
        }


        // Member functions

            virtual scalarField& RobinK()
            {
              return RobinKeff_;
            }

            virtual const scalarField& RobinK() const
            {
              return RobinKeff_;
            }

            virtual scalarField& RobinF()
            {
              return RobinFeff_;
            }

            virtual const scalarField& RobinF() const
            {
              return RobinFeff_;
            }

            //- Evaluate the patch field
            // virtual void evaluate
            // (
            //     const Pstream::commsTypes commsType =
            //         Pstream::commsTypes::blocking
            // );

            virtual void updateCoeffs();

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
