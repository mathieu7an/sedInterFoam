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

#include "nonePPressure.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace granularRheologyModels
{
    defineTypeNameAndDebug(nonePPressure, 0);
    addToRunTimeSelectionTable(PPressureModel, nonePPressure, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::granularRheologyModels::nonePPressure::nonePPressure
(
    const dictionary& dict
)
:
    PPressureModel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::granularRheologyModels::nonePPressure::~nonePPressure()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::granularRheologyModels::nonePPressure::ps
(
    const volScalarField& pf,
    const dimensionedScalar& Bphi,
    const dimensionedScalar& rhos,
    const dimensionedScalar& ds,
    const volScalarField& rhof,
    const volScalarField& nuf,
    const volScalarField& magD,
    const volScalarField& alphas,
    const dimensionedScalar& alphasMax,
    const dimensionedScalar& Alphassmall
) const
{
    // No shear induced pressure
    return scalar(0.0)*pf;
}

Foam::tmp<Foam::volScalarField> Foam::granularRheologyModels::
nonePPressure::alphasEq
(
    const volScalarField& ps,
    const dimensionedScalar& Bphi,
    const dimensionedScalar& rhos,
    const dimensionedScalar& ds,
    const volScalarField& rhof,
    const volScalarField& nuf,
    const volScalarField& magD,
    const dimensionedScalar& alphasMax
) const
{
    // No shear induced pressure
    return scalar(0.0)*magD/magD;
}
// ************************************************************************* //
