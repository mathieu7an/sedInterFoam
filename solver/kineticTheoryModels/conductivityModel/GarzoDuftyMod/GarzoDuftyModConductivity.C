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

#include "GarzoDuftyModConductivity.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(GarzoDuftyModConductivity, 0);

    addToRunTimeSelectionTable
    (
        conductivityModel,
        GarzoDuftyModConductivity,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::GarzoDuftyModConductivity::GarzoDuftyModConductivity
(
const dictionary& dict
)
:
    conductivityModel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::GarzoDuftyModConductivity::~GarzoDuftyModConductivity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::GarzoDuftyModConductivity::kappa
(
    const volScalarField& alphas,
    const volScalarField& Theta,
    const volScalarField& g0,
    const dimensionedScalar& rhos,
    const dimensionedScalar& ds,
    const dimensionedScalar& e
) const
{
    const scalar sqrtPi = sqrt(constant::mathematical::pi);
    const scalar Pi = constant::mathematical::pi;

    return rhos*ds*sqrt(Theta)*225*sqrtPi/1152*
    (
     //Kinetic conductivity
     32*(576./(225*sqrtPi)*alphas + 3./5*pow(1+e, 2)*(2*e-1)*alphas*g0)/
     ((16-7*(1-e))*(1+e)*g0) +
     //Contact conductivity
     32*(576./(225*sqrtPi)*alphas + 3./5*pow(1+e, 2)*(2*e-1)*alphas*g0)/
     ((16-7*(1-e))*(1+e)*g0)* 6./5*(1+e)*alphas*g0 +
     //Bulk conductivity
     2304./(225*Pi)*(1+e)*pow(alphas, 2)*g0
    );
}


// ************************************************************************* //
