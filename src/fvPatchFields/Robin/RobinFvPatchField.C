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

#include "RobinFvPatchField.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::RobinFvPatchField<Type>::RobinFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(p, iF),
    RobinD_(p.size()),
    RobinK_(p.size()),
    RobinF_(p.size())
{
}


template<class Type>
Foam::RobinFvPatchField<Type>::RobinFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchField<Type>(p, iF, dict),
    RobinD_(p.size()),
    RobinK_(p.size()),
    RobinF_(p.size())
{
  if (dict.found("RobinD"))
  {
    RobinD_ = scalarField("RobinD",dict,p.size());
  }
  else
  {
    RobinD_ = scalarField(p.size(),1.0);
  }
  if (dict.found("RobinK"))
  {
    RobinK_ = scalarField("RobinK",dict,p.size());
  }
  else
  {
    RobinK_ = scalarField(p.size(),0.0);
  }
  if (dict.found("RobinF"))
  {
    RobinF_ = Field<Type>("RobinF",dict,p.size());
  }
  else
  {
    RobinF_ = Field<Type>(p.size(),pTraits<Type>::zero);
  }
//    updateCoeffs();
//    evaluate();
}


template<class Type>
Foam::RobinFvPatchField<Type>::RobinFvPatchField
(
    const RobinFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper,
    const bool mappingRequired
)
:
    fvPatchField<Type>(ptf, p, iF, mapper, mappingRequired),
    RobinD_(mapper(ptf.RobinD_)),
    RobinK_(mapper(ptf.RobinK_)),
    RobinF_(mapper(ptf.RobinF_))
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
Foam::RobinFvPatchField<Type>::RobinFvPatchField
(
    const RobinFvPatchField<Type>& ptf
)
:
    fvPatchField<Type>(ptf),
    RobinD_(ptf.RobinD_),
    RobinK_(ptf.RobinK_),
    RobinF_(ptf.RobinF_)
{}


template<class Type>
Foam::RobinFvPatchField<Type>::RobinFvPatchField
(
    const RobinFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(ptf, iF),
    RobinD_(ptf.RobinD_),
    RobinK_(ptf.RobinK_),
    RobinF_(ptf.RobinF_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::RobinFvPatchField<Type>::autoMap
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
void Foam::RobinFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fvPatchField<Type>::rmap(ptf, addr);

    const RobinFvPatchField<Type>& mptf =
        refCast<const RobinFvPatchField<Type>>(ptf);

    RobinD_.rmap(mptf.RobinD_, addr);
    RobinK_.rmap(mptf.RobinK_, addr);
    RobinF_.rmap(mptf.RobinF_, addr);
}

template<class Type>
void Foam::RobinFvPatchField<Type>::evaluate(const Pstream::commsTypes)
{

    if (!this->updated())
    {
        this->updateCoeffs();
    }

    Field<Type>::operator=
    (
        (
          this->patchInternalField() * RobinD() * this->patch().deltaCoeffs()
          + RobinF()
        )
        /
        (
           RobinD() * this->patch().deltaCoeffs()
          - RobinK()
        )
    );

    fvPatchField<Type>::evaluate();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::RobinFvPatchField<Type>::snGrad() const
{
    return
        (
          RobinK() * this->patchInternalField()
          + RobinF()
        ) * this->patch().deltaCoeffs()
        /
        (
          RobinD() * this->patch().deltaCoeffs()
         - RobinK()
        );
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::RobinFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return Type(pTraits<Type>::one) *
    (
        RobinD() * this->patch().deltaCoeffs() )
        /
        (
            RobinD()*this->patch().deltaCoeffs() - RobinK()
        );

}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::RobinFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return
         (
           RobinF()         /
           (
            RobinD() * this->patch().deltaCoeffs()
           - RobinK()
           )
         );
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::RobinFvPatchField<Type>::gradientInternalCoeffs() const
{
    return Type(pTraits<Type>::one) *
    (
      RobinK()*this->patch().deltaCoeffs()
      /
      (
        RobinD() * this->patch().deltaCoeffs()
       - RobinK()
      )
    );
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::RobinFvPatchField<Type>::gradientBoundaryCoeffs() const
{
    return
       (
         RobinF()* this->patch().deltaCoeffs()
         /
         (
          RobinD() * this->patch().deltaCoeffs()
         - RobinK()
         )
       );
}


template<class Type>
void Foam::RobinFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    writeEntry(os,"RobinD",RobinD_);//RobinD_.writeEntry("RobinD", os);
    writeEntry(os,"RobinK",RobinK_);//RobinK_.writeEntry("RobinK", os);
    writeEntry(os,"RobinF",RobinF_);//RobinF_.writeEntry("RobinF", os);
    writeEntry(os,"value",*this);//this->writeEntry("value", os);
}



template<class Type>
void Foam::RobinFvPatchField<Type>::updateCoeffs()
 {

     if (this->updated())
     {
         return;
     }

     fvPatchField<Type>::updateCoeffs();
 }

// ************************************************************************* //
