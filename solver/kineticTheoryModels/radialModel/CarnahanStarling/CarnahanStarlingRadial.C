/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "CarnahanStarlingRadial.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(CarnahanStarlingRadial, 0);

    addToRunTimeSelectionTable
    (
        radialModel,
        CarnahanStarlingRadial,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CarnahanStarlingRadial::CarnahanStarlingRadial(const dictionary& dict)
:
    radialModel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::CarnahanStarlingRadial::~CarnahanStarlingRadial()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::CarnahanStarlingRadial::g0
(
    const volScalarField& alphas,
    const dimensionedScalar& alphasMax
) const
{
    return 1.0/(1.0 - alphas)
         + 3.0*alphas/(2.0*sqr(1.0 - alphas))
         + sqr(alphas)/(2.0*pow(1.0 - alphas, 3));
}


Foam::tmp<Foam::volScalarField> Foam::CarnahanStarlingRadial::g0prime
(
    const volScalarField& alphas,
    const dimensionedScalar& alphasMax
) const
{
    // modified by C.Z.
    return 1.0/sqr(1.0 - alphas)
         + (3.0*(1.0 - alphas) + 6.0*alphas)/(2.0*(1.0 - alphas))
         + (2.0*alphas*(1.0 - alphas) + 3.0*pow(alphas, 2))
          /(2.0*pow(1.0 - alphas, 4));
}


// ************************************************************************* //
