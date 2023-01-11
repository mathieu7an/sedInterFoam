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

#include "GidaspowErgunWenYu.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(GidaspowErgunWenYu, 0);

    addToRunTimeSelectionTable
    (
        dragModel,
        GidaspowErgunWenYu,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::GidaspowErgunWenYu::GidaspowErgunWenYu
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

Foam::GidaspowErgunWenYu::~GidaspowErgunWenYu()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::GidaspowErgunWenYu::K
(
    const volScalarField& Ur
) const
{
    volScalarField alphaf
    (
        max(scalar(1) - alphas_, scalar(1.0e-6))
    );

    //    volScalarField bp = pow(alphaf, -2.65);
    volScalarField bp(pow(alphaf, -phases_.hExp()));
    volScalarField Re
    (
        max(alphaf*Ur*phases_.d()*phases_.sF()/nuf_,
                            scalar(1.0e-9))
    );

    volScalarField Cds
    (
        neg(Re - 1000)*(24.0*(1.0 + 0.15*pow(Re, 0.687))/Re)
      + pos(Re - 1000)*0.44
    );

/*
    // Wen and Yu (1966)
    // modified in May 17, 2012, by C.Z.
    //tmp<volScalarField> tKWenYu = (0.75*Cds*phasew_.rho()*Ur*bp
                                   /(phases_.d()*phases_.sF()));
    volScalarField tKWenYu = (0.75*Cds*phasew_.rho()*Ur*bp
                            /(phases_.d()*phases_.sF()));
    volScalarField& KWenYu = tKWenYu();

    // Ergun
    forAll(alphaf, cellj)
    {
        if (alphaf[cellj] < 0.8)
        {
            tKWenYu[cellj] =
                150.0*alphas_[cellj]*phasew_.nu().value()*phasew_.rho().value()
                  /sqr(alphaf[cellj]*phases_.d().value())
                 + 1.75*phasew_.rho().value()*Ur[cellj]
                  /(alphaf[cellj]*phases_.d().value());
        }
    }
// WARNING: remove this line will makes the parallel computations "instable"
    tKWenYu.correctBoundaryConditions();

    return KWenYu;
*/
//        pos0(alphaf - 0.8)*(0.75*Cds*phasew_.rho()*Ur*bp
//      /(phases_.d()*phases_.sF()))    //line for sedfoam-5.0
    return pos(alphaf - 0.8)
          *(0.75*Cds*rhof_*Ur*bp/(phases_.d()*phases_.sF()))
         + neg(alphaf - 0.8)
          *(150.0*alphas_*nuf_*rhof_/sqr(alphaf*phases_.d())
                + 1.75*rhof_*Ur/(alphaf*phases_.d()));
}


// ************************************************************************* //
