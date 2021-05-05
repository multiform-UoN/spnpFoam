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
    Foam::mappedChemicalKineticsSolidFvPatchScalarField

Description
    This boundary condition imposes continuity of flux alongside a given reaction
    rate set to the fluxes.
Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    \endverbatim


SourceFiles
    mappedChemicalKineticsSolidFvPatchScalarField.C

Authors:
    Robert Barnett, Nottingham (2020)
\*---------------------------------------------------------------------------*/

#ifndef mappedChemicalKineticsSolidFvPatchScalarField_H
#define mappedChemicalKineticsSolidFvPatchScalarField_H

//#include "mappedContinuitySolidFvPatchScalarField.H"
#include "RobinFvPatchFields.H"
#include "newMappedPatchFieldBase.H"

#include "mappedChemicalKineticsSolidFvPatchScalarField.H"
#include "reactionRate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    /*---------------------------------------------------------------------------*\
             Class mappedChemicalKineticsSolidFvPatchScalarField Declaration
    \*---------------------------------------------------------------------------*/
    class mappedChemicalKineticsSolidFvPatchScalarField
    :
    public RobinFvPatchScalarField,
    public newMappedPatchFieldBase<scalar>
    {
        private:

            // Local (this patch) Flux variable
            word localPhi_;

            // Mapped (other mapped patch) Flux variable
            word mappedPhi_;

            // Effective RobinK
            scalarField RobinKeff_;

            // Effective RobinF
            scalarField RobinFeff_;

            // Effective Dfluid
            scalarField Dfluid_;

            // Pointer to reaction rate
            autoPtr<reactionRates::reactionRate> R_;

        protected:

        // Protected Member Functions

            const mappedPatchBase& mapper
            (
                const fvPatch& p,
                const DimensionedField<scalar, volMesh>& iF
            );

        public:

            //- Runtime type information
            TypeName("mappedChemicalKineticsSolid");


        // Constructors

            //- Construct from patch and internal field
            mappedChemicalKineticsSolidFvPatchScalarField
            (
                const fvPatch&,
                const DimensionedField<scalar, volMesh>&
            );

            //- Construct from patch, internal field and dictionary
            mappedChemicalKineticsSolidFvPatchScalarField
            (
                const fvPatch&,
                const DimensionedField<scalar, volMesh>&,
                const dictionary&
            );

            //- Construct by mapping given mappedChemicalKineticsSolidFvPatchScalarField
            //  onto a new patch
            mappedChemicalKineticsSolidFvPatchScalarField
            (
                const mappedChemicalKineticsSolidFvPatchScalarField&,
                const fvPatch&,
                const DimensionedField<scalar, volMesh>&,
                const fvPatchFieldMapper&
            );

            //- Construct as copy
            mappedChemicalKineticsSolidFvPatchScalarField
            (
                const mappedChemicalKineticsSolidFvPatchScalarField&
            );

            //- Construct and return a clone
            tmp<fvPatchScalarField> clone() const
            {
                return tmp<fvPatchScalarField>
                (
                    new mappedChemicalKineticsSolidFvPatchScalarField(*this)
                );
            }

            //- Construct as copy setting internal field reference
            mappedChemicalKineticsSolidFvPatchScalarField
            (
                const mappedChemicalKineticsSolidFvPatchScalarField&,
                const DimensionedField<scalar, volMesh>&
            );

            //- Construct and return a clone setting internal field reference
            tmp<fvPatchScalarField> clone
            (
                const DimensionedField<scalar, volMesh>& iF
            ) const
            {
                return tmp<fvPatchScalarField>
                (
                    new mappedChemicalKineticsSolidFvPatchScalarField(*this, iF)
                );
            }

        // Member functions

            //- Evaluate the patch field
            // virtual void evaluate
            // (
            //     const Pstream::commsTypes commsType =
            //         Pstream::commsTypes::blocking
            // );

            //- Map (and resize as needed) from self given a mapping object
            void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            void rmap
            (
                const fvPatchField<scalar>&,
                const labelList&
            );

            virtual scalarField& RobinF()
            {
                return RobinFeff_;
            }

            virtual const scalarField& RobinF() const
            {
                return RobinFeff_;
            }

            virtual scalarField& RobinK()
            {
                return RobinKeff_;
            }

            virtual const scalarField& RobinK() const
            {
                return RobinKeff_;
            }

            void updateCoeffs();

            //- Write
            void write(Ostream&) const;


    };
}



#endif
