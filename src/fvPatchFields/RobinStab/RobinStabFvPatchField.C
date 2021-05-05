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

#include "RobinStabFvPatchField.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::RobinStabFvPatchField<Type>::RobinStabFvPatchField
(
  const fvPatch& p,
  const DimensionedField<Type, volMesh>& iF
)
:
fvPatchField<Type>(p, iF),
RobinD_(p.size()),
RobinK_(p.size()),
RobinF_(p.size()),
beta_(1.)
{
}


template<class Type>
Foam::RobinStabFvPatchField<Type>::RobinStabFvPatchField
(
  const fvPatch& p,
  const DimensionedField<Type, volMesh>& iF,
  const dictionary& dict
)
:
fvPatchField<Type>(p, iF, dict, false),
RobinD_("RobinD", dict, p.size()),
RobinK_("RobinK", dict, p.size()),
RobinF_("RobinF", dict, p.size()),
beta_(readScalar(dict.lookup("beta")))
{
  evaluate();
}


template<class Type>
Foam::RobinStabFvPatchField<Type>::RobinStabFvPatchField
(
  const RobinStabFvPatchField<Type>& ptf,
  const fvPatch& p,
  const DimensionedField<Type, volMesh>& iF,
  const fvPatchFieldMapper& mapper,
  const bool mappingRequired
)
:
fvPatchField<Type>(ptf, p, iF, mapper, mappingRequired),
RobinD_(mapper(ptf.RobinD_)),
RobinK_(mapper(ptf.RobinK_)),
RobinF_(mapper(ptf.RobinF_)),
beta_(ptf.beta_)
{
  if (mappingRequired && notNull(iF) && mapper.hasUnmapped())
  {
    WarningInFunction
    << "On field " << iF.name() << " patch " << p.name()
    << " patchField " << this->type()
    << " : mapper does not map all values." << nl
    << "    To avoid this warning fully specify the mapping in derived"
    << " patch fields." << endl;
  }
}


template<class Type>
Foam::RobinStabFvPatchField<Type>::RobinStabFvPatchField
(
  const RobinStabFvPatchField<Type>& ptf
)
:
fvPatchField<Type>(ptf),
RobinD_(ptf.RobinD_),
RobinK_(ptf.RobinK_),
RobinF_(ptf.RobinF_),
beta_(ptf.beta_)
{}


  template<class Type>
  Foam::RobinStabFvPatchField<Type>::RobinStabFvPatchField
  (
    const RobinStabFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
  )
  :
  fvPatchField<Type>(ptf, iF),
  RobinD_(ptf.RobinD_),
  RobinK_(ptf.RobinK_),
  RobinF_(ptf.RobinF_),
  beta_(ptf.beta_)
  {}


    // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

    template<class Type>
    void Foam::RobinStabFvPatchField<Type>::autoMap
    (
      const fvPatchFieldMapper& m
    )
    {
      fvPatchField<Type>::autoMap(m);
      m(RobinD_,RobinD_);//.autoMap(m);
      m(RobinK_,RobinK_);//RobinK_.autoMap(m);
      m(RobinF_,RobinF_);//.autoMap(m);
    }


    template<class Type>
    void Foam::RobinStabFvPatchField<Type>::rmap
    (
      const fvPatchField<Type>& ptf,
      const labelList& addr
    )
    {
      fvPatchField<Type>::rmap(ptf, addr);

      const RobinStabFvPatchField<Type>& mptf =
      refCast<const RobinStabFvPatchField<Type>>(ptf);

      RobinD_.rmap(mptf.RobinD_, addr);
      RobinK_.rmap(mptf.RobinK_, addr);
      RobinF_.rmap(mptf.RobinF_, addr);
    }

    template<class Type>
    void Foam::RobinStabFvPatchField<Type>::evaluate(const Pstream::commsTypes)
    {
      if (!this->updated())
      {
        this->updateCoeffs();
      }

      Field<Type>::operator=
      (
        (
          this->patchInternalField() * RobinD() * this->patch().deltaCoeffs()
          + RobinFStab()
        )
        /
        (
          RobinD() * this->patch().deltaCoeffs()
          + RobinKStab()
        )
      );

      fvPatchField<Type>::evaluate();
    }


    template<class Type>
    Foam::tmp<Foam::Field<Type>>
    Foam::RobinStabFvPatchField<Type>::snGrad() const
    {
      return
      (
        (
          RobinKStab() * this->patchInternalField()
          +
          RobinFStab()
        ) * this->patch().deltaCoeffs()
        /
        (
          RobinD() * this->patch().deltaCoeffs()
          - RobinKStab()
        )
      );
    }


    template<class Type>
    Foam::tmp<Foam::Field<Type>>
    Foam::RobinStabFvPatchField<Type>::valueInternalCoeffs
    (
      const tmp<scalarField>&
    ) const
    {
      return Type(pTraits<Type>::one) *
      (
        RobinD() * this->patch().deltaCoeffs() )
        /
        (
          RobinD()*this->patch().deltaCoeffs() - RobinKStab()
        );

      }


      template<class Type>
      Foam::tmp<Foam::Field<Type>>
      Foam::RobinStabFvPatchField<Type>::valueBoundaryCoeffs
      (
        const tmp<scalarField>&
      ) const
      {
        return
        (
          RobinFStab()
          /
          (
            RobinD() * this->patch().deltaCoeffs()
            - RobinKStab()
          )
        );
      }


      template<class Type>
      Foam::tmp<Foam::Field<Type>>
      Foam::RobinStabFvPatchField<Type>::gradientInternalCoeffs() const
      {
        return Type(pTraits<Type>::one) *
        (
          RobinK()*this->patch().deltaCoeffs()
          /
          (
            RobinD() * this->patch().deltaCoeffs()
            - RobinKStab()
          )
        );
      }


      template<class Type>
      Foam::tmp<Foam::Field<Type>>
      Foam::RobinStabFvPatchField<Type>::gradientBoundaryCoeffs() const
      {
        return
        (
          (
            RobinFStab()
            *
            this->patch().deltaCoeffs()
          )
          /
          (
            RobinD() * this->patch().deltaCoeffs()
            - RobinKStab()
          )
        );
      }


      template<class Type>
      void Foam::RobinStabFvPatchField<Type>::write(Ostream& os) const
      {
        fvPatchField<Type>::write(os);
        writeEntry(os,"RobinD",RobinD_);//RobinD_.writeEntry("RobinD", os);
        writeEntry(os,"RobinK",RobinK_);//RobinK_.writeEntry("RobinK", os);
        writeEntry(os,"RobinF",RobinF_);//RobinF_.writeEntry("RobinF", os);
        scalarField KStab(RobinKStab());
        writeEntry(os,"RobinKStab",KStab);
        writeEntry(os,"beta",beta());
        writeEntry(os,"value",*this);//this->writeEntry("value", os);
      }



      template<class Type>
      void Foam::RobinStabFvPatchField<Type>::updateCoeffs()
      {
        if (this->updated())
        {
          return;
        }

        fvPatchField<Type>::updateCoeffs();

      }

      // ************************************************************************* //
