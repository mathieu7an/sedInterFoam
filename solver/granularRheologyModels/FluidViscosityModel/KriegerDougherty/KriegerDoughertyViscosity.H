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
    Foam::KriegerDoughertyViscosity

Description

    This file is in the \dir FluidViscosityModel/KriegerDougherty subdirectory
    inside the \dir FluidViscosityModel subdirectory

SourceFiles
    KriegerDoughertyViscosity.C

\*---------------------------------------------------------------------------*/

#ifndef KriegerDoughertyViscosity_H
#define KriegerDoughertyViscosity_H

#include "FluidViscosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace granularRheologyModels
{

/*---------------------------------------------------------------------------*\
                           Class KriegerDoughertyViscosity Declaration
\*---------------------------------------------------------------------------*/

class KriegerDoughertyViscosity
:
    public FluidViscosityModel
{

public:

    //- Runtime type information
    TypeName("KriegerDougherty");


    // Constructors

    //- Construct from components
    KriegerDoughertyViscosity(const dictionary& dict);


    //- Destructor
    virtual ~KriegerDoughertyViscosity();


    // Member functions

    tmp<volScalarField> nuvf
    (
        const volScalarField& alphas,
        const volScalarField& nuf,
        const dimensionedScalar& alphasMax,
        const dimensionedScalar& Alphassmall,
        const dimensionedScalar& n
    ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace granularRheologyModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
