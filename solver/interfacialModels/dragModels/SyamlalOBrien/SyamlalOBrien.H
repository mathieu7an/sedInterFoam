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
    Foam::SyamlalOBrien

Description
    Syamlal, M., Rogers, W. and O'Brien, T. J. (1993) MFIX documentation,
    Theory Guide. Technical Note DOE/METC-94/1004. Morgantown, West Virginia,
    USA.

SourceFiles
    SyamlalOBrien.C

\*---------------------------------------------------------------------------*/

#ifndef SyamlalOBrien_H
#define SyamlalOBrien_H

#include "dragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class SyamlalOBrien Declaration
\*---------------------------------------------------------------------------*/

class SyamlalOBrien
:
    public dragModel
{

public:

    //- Runtime type information
    TypeName("SyamlalOBrien");


    // Constructors

        //- Construct from components
        SyamlalOBrien
        (
            const dictionary& interfaceDict,
            const phaseModel& phases,
            const volScalarField& rhof,
            const volScalarField& nuf
        );


    //- Destructor
    virtual ~SyamlalOBrien();


    // Member Functions

        tmp<volScalarField> K(const volScalarField& Ur) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
