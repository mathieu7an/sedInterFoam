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

#include "DallaValle.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(DallaValle, 0);

    addToRunTimeSelectionTable
    (
        dragModel,
        DallaValle,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::DallaValle::DallaValle
(
    const dictionary& interfaceDict,
    const phaseModel& phases,
    const volScalarField& rhof,
    const volScalarField& nuf
)
:
    dragModel(interfaceDict, phases, rhof, nuf)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::DallaValle::~DallaValle()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::DallaValle::K
(
    const volScalarField& Ur
) const
{
    volScalarField alphaf(max(scalar(1) - alphas_, scalar(1e-6)));
    volScalarField bp(pow(alphaf, -phases_.hExp()-1));
    volScalarField Re
    (
        max(Ur*phases_.d()*phases_.sF()/nuf_, scalar(1.0e-9))
    );

    volScalarField Cds
    (
      0.4+24.4/Re
    );

    return 0.75*Cds*rhof_*Ur*bp/(phases_.d()*phases_.sF());
}


// ************************************************************************* //
