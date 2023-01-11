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
\*---------------------------------------------------------------------------*/

#include "EinsteinViscosity.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace granularRheologyModels
    {
        defineTypeNameAndDebug(EinsteinViscosity, 0);
        addToRunTimeSelectionTable
        (
            FluidViscosityModel, EinsteinViscosity, dictionary
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::granularRheologyModels::EinsteinViscosity::EinsteinViscosity
(
    const dictionary& dict
):
    FluidViscosityModel(dict)
{}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::granularRheologyModels::EinsteinViscosity::~EinsteinViscosity()
{}
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::granularRheologyModels::
                                      EinsteinViscosity::nuvf
(
    const volScalarField& alphas,
    const volScalarField& nuf,
    const dimensionedScalar& alphasMax,
    const dimensionedScalar& Alphassmall,
    const dimensionedScalar& n
) const
{
    return nuf*(1.0 + 2.5 * alphas)/(1-alphas);
}


// ************************************************************************* //
