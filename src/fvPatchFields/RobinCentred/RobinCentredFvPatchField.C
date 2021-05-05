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

#include "RobinCentredFvPatchField.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::RobinCentredFvPatchField<Type>::RobinCentredFvPatchField
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
Foam::RobinCentredFvPatchField<Type>::RobinCentredFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchField<Type>(p, iF, dict, false),
    RobinD_("RobinD", dict, p.size()),
    RobinK_("RobinK", dict, p.size()),
    RobinF_("RobinF", dict, p.size())
{
    evaluate();
}


template<class Type>
Foam::RobinCentredFvPatchField<Type>::RobinCentredFvPatchField
(
    const RobinCentredFvPatchField<Type>& ptf,
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
Foam::RobinCentredFvPatchField<Type>::RobinCentredFvPatchField
(
    const RobinCentredFvPatchField<Type>& ptf
)
:
    fvPatchField<Type>(ptf),
    RobinD_(ptf.RobinD_),
    RobinK_(ptf.RobinK_),
    RobinF_(ptf.RobinF_)
{}


template<class Type>
Foam::RobinCentredFvPatchField<Type>::RobinCentredFvPatchField
(
    const RobinCentredFvPatchField<Type>& ptf,
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
void Foam::RobinCentredFvPatchField<Type>::autoMap
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
void Foam::RobinCentredFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fvPatchField<Type>::rmap(ptf, addr);

    const RobinCentredFvPatchField<Type>& mptf =
        refCast<const RobinCentredFvPatchField<Type>>(ptf);

    RobinD_.rmap(mptf.RobinD_, addr);
    RobinK_.rmap(mptf.RobinK_, addr);
    RobinF_.rmap(mptf.RobinF_, addr);
}

template<class Type>
void Foam::RobinCentredFvPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    Field<Type>::operator=
    (
        (
          this->patchInternalField() *
          ( RobinK() + this->patch().deltaCoeffs()*RobinD() )
          + RobinF()
        )
        /
        (
           RobinD() * this->patch().deltaCoeffs()
        )
    );

    fvPatchField<Type>::evaluate();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::RobinCentredFvPatchField<Type>::snGrad() const
{
    return
    (
      this->patchInternalField() * RobinK()
      + RobinF()
    )
    /
    (
       RobinD()
    );
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::RobinCentredFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return Type(pTraits<Type>::one) *
      ( RobinK() + this->patch().deltaCoeffs()*RobinD() )
      /
      (
        RobinD() * this->patch().deltaCoeffs()
      );
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::RobinCentredFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return
         (
           RobinF()
           /
           (
            RobinD() * this->patch().deltaCoeffs()
           )
         );
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::RobinCentredFvPatchField<Type>::gradientInternalCoeffs() const
{
    return Type(pTraits<Type>::one) *
    (
      RobinK()
      /
      RobinD()
    );
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::RobinCentredFvPatchField<Type>::gradientBoundaryCoeffs() const
{
    return
       (
         RobinF()
         /
         RobinD()
       );
}


template<class Type>
void Foam::RobinCentredFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    writeEntry(os,"RobinD",RobinD_);//RobinD_.writeEntry("RobinD", os);
    writeEntry(os,"RobinK",RobinK_);//RobinK_.writeEntry("RobinK", os);
    writeEntry(os,"RobinF",RobinF_);//RobinF_.writeEntry("RobinF", os);
    writeEntry(os,"value",*this);//this->writeEntry("value", os);
}



template<class Type>
void Foam::RobinCentredFvPatchField<Type>::updateCoeffs()
 {
     if (this->updated())
     {
         return;
     }

     fvPatchField<Type>::updateCoeffs();
 }

// ************************************************************************* //
