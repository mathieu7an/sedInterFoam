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

#include "TorquatoRadial.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(TorquatoRadial, 0);

    addToRunTimeSelectionTable
    (
        radialModel,
        TorquatoRadial,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::TorquatoRadial::TorquatoRadial(const dictionary& dict)
:
    radialModel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::TorquatoRadial::~TorquatoRadial()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::TorquatoRadial::g0
(
    const volScalarField& alphas,
    const dimensionedScalar& alphasMax
) const
{
    return (neg(alphas-0.49)*(2-alphas)/(2*pow(1-alphas, 3)) +
     pos(alphas-0.49)*(2-0.49)*(alphasMax-0.49)/
     (2*pow(1-0.49, 3)*(alphasMax-alphas)));
}


Foam::tmp<Foam::volScalarField> Foam::TorquatoRadial::g0prime
(
    const volScalarField& alphas,
    const dimensionedScalar& alphasMax
) const
{
    // need to be updated: actual version is CarnahanStarling
    // Never used in the actual code
    return 1.0/sqr(1.0 - alphas)
         + (3.0*(1.0 - alphas) + 6.0*alphas)/(2.0*(1.0 - alphas))
         + (2.0*alphas*(1.0 - alphas) + 3.0*pow(alphas, 2))
          /(2.0*pow(1.0 - alphas, 4));
}


// ************************************************************************* //
