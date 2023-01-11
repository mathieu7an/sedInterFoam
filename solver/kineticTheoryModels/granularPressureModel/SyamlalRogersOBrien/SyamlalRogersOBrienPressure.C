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

#include "SyamlalRogersOBrienPressure.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(SyamlalRogersOBrienPressure, 0);

    addToRunTimeSelectionTable
    (
        granularPressureModel,
        SyamlalRogersOBrienPressure,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SyamlalRogersOBrienPressure::SyamlalRogersOBrienPressure
(
    const dictionary& dict
)
:
    granularPressureModel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::SyamlalRogersOBrienPressure::~SyamlalRogersOBrienPressure()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::SyamlalRogersOBrienPressure::
granularPressureCoeff
(
    const volScalarField& alphas,
    const volScalarField& g0,
    const dimensionedScalar& rhos,
    const dimensionedScalar& e
) const
{

    return 2.0*rhos*(1.0 + e)*sqr(alphas)*g0;
}


Foam::tmp<Foam::volScalarField> Foam::SyamlalRogersOBrienPressure::
granularPressureCoeffPrime
(
    const volScalarField& alphas,
    const volScalarField& g0,
    const volScalarField& g0prime,
    const dimensionedScalar& rhos,
    const dimensionedScalar& e
) const
{
    return rhos*alphas*(1.0 + e)*(4.0*g0 + 2.0*g0prime*alphas);
}


// ************************************************************************* //
