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

Class
    Foam::dragModel

Description

SourceFiles
    dragModel.C
    newDragModel.C

\*---------------------------------------------------------------------------*/

#ifndef dragModel_H
#define dragModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dictionary.H"
#include "phaseModel.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dragModel Declaration
\*---------------------------------------------------------------------------*/

class dragModel
{
protected:

    // Protected data

    const dictionary& interfaceDict_;
    const volScalarField& alphas_;
    const phaseModel& phases_;
    const volScalarField& rhof_;
    const volScalarField& nuf_;


public:

    //- Runtime type information
    TypeName("dragModel");


    // Declare runtime construction

    declareRunTimeSelectionTable
    (
        autoPtr,
        dragModel,
        dictionary,
        (
            const dictionary& interfaceDict,
            const phaseModel& phases,
            const volScalarField& rhof,
            const volScalarField& nuf
        ),
        (interfaceDict, phases, rhof, nuf)
    );


    // Constructors

    dragModel
    (
        const dictionary& interfaceDict,
        const phaseModel& phases,
        const volScalarField& rhof,
        const volScalarField& nuf
    );


    //- Destructor
    virtual ~dragModel();


    // Selectors

    static autoPtr<dragModel> New
    (
        const dictionary& interfaceDict,
        const phaseModel& phases,
        const volScalarField& rhof,
        const volScalarField& nuf
    );


    // Member Functions

    //- the dragfunction K used in the momentum eq.
    //    ddt(alphas.rhos.Us) + ... = ... alphas.alphaw.K.(Us-Uw)
    //    ddt(alphaw.rhow.Uw) + ...  = ... alphas.alphaw.K.(Uw-Us)
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
    // for numerical reasons alphas and alphaw has been
    // extracted from the dragFunction K,
    // so you MUST divide K by alphas.alphaflu when implemnting the drag function
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
    virtual tmp<volScalarField> K(const volScalarField& Ur) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
