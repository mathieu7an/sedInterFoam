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

#include "LunPressure.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(LunPressure, 0);

    addToRunTimeSelectionTable
    (
        granularPressureModel,
        LunPressure,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LunPressure::LunPressure(const dictionary& dict)
:
    granularPressureModel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::LunPressure::~LunPressure()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::LunPressure::granularPressureCoeff
(
    const volScalarField& alphas,
    const volScalarField& g0,
    const dimensionedScalar& rhos,
    const dimensionedScalar& e
) const
{

    return rhos*alphas*(1.0 + 2.0*(1.0 + e)*alphas*g0);
}


Foam::tmp<Foam::volScalarField> Foam::LunPressure::granularPressureCoeffPrime
(
    const volScalarField& alphas,
    const volScalarField& g0,
    const volScalarField& g0prime,
    const dimensionedScalar& rhos,
    const dimensionedScalar& e
) const
{
    return rhos*(1.0 + alphas*(1.0 + e)*(4.0*g0 + 2.0*g0prime*alphas));
}

// ************************************************************************* //
