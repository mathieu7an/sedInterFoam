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

#include "ChialvoSundaresanFricRadial.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ChialvoSundaresanFricRadial, 0);

    addToRunTimeSelectionTable
    (
        radialModel,
        ChialvoSundaresanFricRadial,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ChialvoSundaresanFricRadial::ChialvoSundaresanFricRadial
(
 const dictionary& dict
)
:
    radialModel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ChialvoSundaresanFricRadial::~ChialvoSundaresanFricRadial()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::ChialvoSundaresanFricRadial::g0
(
    const volScalarField& alphas,
    const dimensionedScalar& alphasMax
) const
{
    return (2-alphas)/(2*pow((1-alphas), 3)) +
     2.71*pow(alphas, 2)/pow(alphasMax-alphas, 1.5);
}


Foam::tmp<Foam::volScalarField> Foam::ChialvoSundaresanFricRadial::g0prime
(
    const volScalarField& alphas,
    const dimensionedScalar& alphasMax
) const
{
    return 3*1.7/(2*alphasMax)*pow(mag(1.0 - alphas/alphasMax), -2.7);
}


// ************************************************************************* //
