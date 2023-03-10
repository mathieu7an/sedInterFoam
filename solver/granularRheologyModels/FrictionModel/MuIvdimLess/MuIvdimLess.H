/*---------------------------------------------------------------------------*\
Copyright (C) 2015 Cyrille Bonamy, Julien Chauchat, Tian-Jian Hsu
                   and contributors

License
    This file is part of SedFOAM.

    SedFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SedFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with SedFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::MuIvdimLess

Description
    MuI rheology in dimensionless form. Lengths are made dimensionless by d the
    particle diameter and velocity by sqrt(g d).

    This file is in the \dir FrictionModel/MuIvdimLess subdirectory
    inside the \dir FrictionModel subdirectory

SourceFiles
    MuIvdimLess.C

\*---------------------------------------------------------------------------*/

#ifndef MuIvdimLess_H
#define MuIvdimLess_H

#include "FrictionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace granularRheologyModels
{

/*---------------------------------------------------------------------------*\
                           Class MuIvdimLess Declaration
\*---------------------------------------------------------------------------*/

class MuIvdimLess
:
    public FrictionModel
{

public:

    //- Runtime type information
    TypeName("MuIvdimLess");


    // Constructors

    //- Construct from components
    MuIvdimLess(const dictionary& dict);


    //- Destructor
    virtual ~MuIvdimLess();


    // Member functions

    tmp<volScalarField> muI
    (
        const dimensionedScalar& muss,
        const dimensionedScalar& mu2,
        const dimensionedScalar& I0,
        const volScalarField& ps,
        const dimensionedScalar& rhos,
        const dimensionedScalar& ds,
        const volScalarField& rhof,
        const volScalarField& nuf,
        const volScalarField& magD,
        const dimensionedScalar& Dsmall
    ) const;


    tmp<volScalarField> I
    (
        const volScalarField& ps,
        const dimensionedScalar& rhos,
        const dimensionedScalar& ds,
        const volScalarField& rhof,
        const volScalarField& nuf,
        const volScalarField& magD
    ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace granularRheologyModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
