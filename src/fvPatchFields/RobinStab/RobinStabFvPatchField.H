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
Foam::RobinStabFvPatchField

Group
grpGenericBoundaryConditions

Description
This boundary condition provides a base class for 'RobinStab' type boundary
conditions
\f[
D \nabla_\perp x  = K x + F
\f]

where
\vartable
x_p   | patch values
x_c   | patch internal cell values
K     | reaction term
\Delta| inverse distance from face centre to internal cell centre
F     | Forcing term
D     | Diffusion
\endvartable


Usage
\table
Property     | Description             | Required    | Default value
RobinD       | Diffusion               | yes         |
RobinF       | Source Term             | yes         |
RobinStabR       | Reaction                | yes         |
\endtable

SourceFiles
RobinStabFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef RobinStabFvPatchField_H
#define RobinStabFvPatchField_H

#include "fvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

  /*---------------------------------------------------------------------------*\
  Class RobinStabFvPatchField Declaration
  \*---------------------------------------------------------------------------*/

  template<class Type>
  class RobinStabFvPatchField
  :
  public fvPatchField<Type>
  {
    // Private data

    //- Diffusion field
    scalarField RobinD_;

    //- Reaction field
    scalarField RobinK_;

    //- Source term
    Field<Type> RobinF_;

    //-Stabilisation parameter
    scalar beta_;

  public:

    //- Runtime type information
    TypeName("RobinStab");


    // Constructors

    //- Construct from patch and internal field
    RobinStabFvPatchField
    (
      const fvPatch&,
      const DimensionedField<Type, volMesh>&
    );

    //- Construct from patch, internal field and dictionary
    RobinStabFvPatchField
    (
      const fvPatch&,
      const DimensionedField<Type, volMesh>&,
      const dictionary&
    );

    //- Construct by mapping the given RobinStabFvPatchField onto a new patch
    RobinStabFvPatchField
    (
      const RobinStabFvPatchField<Type>&,
      const fvPatch&,
      const DimensionedField<Type, volMesh>&,
      const fvPatchFieldMapper&,
      const bool mappingRequired=true
    );

    //- Construct as copy
    RobinStabFvPatchField
    (
      const RobinStabFvPatchField<Type>&
    );

    //- Construct and return a clone
    virtual tmp<fvPatchField<Type>> clone() const
    {
      return tmp<fvPatchField<Type>>
      (
        new RobinStabFvPatchField<Type>(*this)
      );
    }

    //- Construct as copy setting internal field reference
    RobinStabFvPatchField
    (
      const RobinStabFvPatchField<Type>&,
      const DimensionedField<Type, volMesh>&
    );

    //- Construct and return a clone setting internal field reference
    virtual tmp<fvPatchField<Type>> clone
    (
      const DimensionedField<Type, volMesh>& iF
    ) const
    {
      return tmp<fvPatchField<Type>>
      (
        new RobinStabFvPatchField<Type>(*this, iF)
      );
    }


    // Member functions

    // Access

    //- Return true if this patch field fixes a value.
    //  Needed to check if a level has to be specified while solving
    //  Poissons equations.
    virtual bool fixesValue() const
    {
      return false;
    }

    //- Return false: this patch field is not altered by assignment
    virtual bool assignable() const
    {
      return true;
    }


    // Return defining fields

    virtual scalarField& RobinD()
    {
      return RobinD_;
    }

    virtual const scalarField& RobinD() const
    {
      return RobinD_;
    }

    virtual scalarField& RobinK()
    {
      return RobinK_;
    }

    virtual const scalarField& RobinK() const
    {
      return RobinK_;
    }

    virtual Field<Type>& RobinF()
    {
      return RobinF_;
    }

    virtual const Field<Type>& RobinF() const
    {
      return RobinF_;
    }

    virtual const scalar& beta() const
    {
      return beta_;
    }

    virtual scalar& beta()
    {
      return beta_;
    }

    //- Return stabilised coefficients
    tmp<scalarField> RobinKStab()
    {
      //        return -min(beta()*RobinD(),mag(RobinK()));
      return min(beta()*RobinD(),RobinK());
    }

    tmp<Field<Type>> RobinFStab()
    {
      //        return RobinF() + (RobinKStab()-RobinK())*(*this);
      return RobinF() + (-RobinKStab()+RobinK())*(*this);
    }

    const tmp<scalarField> RobinKStab() const
    {
      //        return -min(beta()*RobinD(),mag(RobinK()));
      return min(beta()*RobinD(),RobinK());
    }

    const tmp<Field<Type>> RobinFStab() const
    {
      //        return RobinF() + (RobinKStab()-RobinK())*(*this);
      return RobinF() + (-RobinKStab()+RobinK())*(*this);
    }

    // Mapping functions

    //- Map (and resize as needed) from self given a mapping object
    virtual void autoMap
    (
      const fvPatchFieldMapper&
    );

    //- Reverse map the given fvPatchField onto this fvPatchField
    virtual void rmap
    (
      const fvPatchField<Type>&,
      const labelList&
    );


    // Evaluation functions

    //- Return gradient at boundary
    virtual tmp<Field<Type>> snGrad() const;

    //- Evaluate the patch field
    virtual void evaluate
    (
      const Pstream::commsTypes commsType =
      Pstream::commsTypes::blocking
    );

    //- Return the matrix diagonal coefficients corresponding to the
    //  evaluation of the value of this patchField with given weights
    virtual tmp<Field<Type>> valueInternalCoeffs
    (
      const tmp<scalarField>&
    ) const;

    //- Return the matrix source coefficients corresponding to the
    //  evaluation of the value of this patchField with given weights
    virtual tmp<Field<Type>> valueBoundaryCoeffs
    (
      const tmp<scalarField>&
    ) const;

    //- Return the matrix diagonal coefficients corresponding to the
    //  evaluation of the gradient of this patchField
    virtual tmp<Field<Type>> gradientInternalCoeffs() const;

    //- Return the matrix source coefficients corresponding to the
    //  evaluation of the gradient of this patchField
    virtual tmp<Field<Type>> gradientBoundaryCoeffs() const;

    virtual void updateCoeffs();


    //- Write
    virtual void write(Ostream&) const;


    // Member operators

    virtual void operator=(const UList<Type>&) {}

    virtual void operator=(const fvPatchField<Type>&) {}
    virtual void operator+=(const fvPatchField<Type>&) {}
    virtual void operator-=(const fvPatchField<Type>&) {}
    virtual void operator*=(const fvPatchField<scalar>&) {}
    virtual void operator/=(const fvPatchField<scalar>&) {}

    virtual void operator+=(const Field<Type>&) {}
    virtual void operator-=(const Field<Type>&) {}

    virtual void operator*=(const Field<scalar>&) {}
    virtual void operator/=(const Field<scalar>&) {}

    virtual void operator=(const Type&) {}
    virtual void operator+=(const Type&) {}
    virtual void operator-=(const Type&) {}
    virtual void operator*=(const scalar) {}
    virtual void operator/=(const scalar) {}
  };


  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "RobinStabFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
