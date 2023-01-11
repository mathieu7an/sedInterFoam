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

#include "MuIvdimLess.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace granularRheologyModels
{
    defineTypeNameAndDebug(MuIvdimLess, 0);
    addToRunTimeSelectionTable(FrictionModel, MuIvdimLess, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::granularRheologyModels::MuIvdimLess::MuIvdimLess(const dictionary& dict)
:
    FrictionModel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::granularRheologyModels::MuIvdimLess::~MuIvdimLess()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::granularRheologyModels::MuIvdimLess::muI
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
) const
{
    return muss + (mu2 - muss)*magD/(I0*rhos/rhof*ds*ps/nuf + magD + Dsmall);
}


Foam::tmp<Foam::volScalarField> Foam::granularRheologyModels::MuIvdimLess::I
(
    const volScalarField& ps,
    const dimensionedScalar& rhos,
    const dimensionedScalar& ds,
    const volScalarField& rhof,
    const volScalarField& nuf,
    const volScalarField& magD
) const
{
    return rhof*nuf*magD / (rhos*ds*ps);
}

// ************************************************************************* //
