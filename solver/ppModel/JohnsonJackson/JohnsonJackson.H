/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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
    Foam::JohnsonJackson

Description
    Johnson, P., & Jackson, R. (1987).
    Frictional–collisional constitutive relations for granular materials,
    with application to plane shearing.
    Journal of Fluid Mechanics, 176, 67-93. doi:10.1017/S0022112087000570

This file is in the \dir ppModel/JohnsonJackson directory.

SourceFiles
    JohnsonJackson.C

\*---------------------------------------------------------------------------*/

#ifndef JohnsonJackson_H
#define JohnsonJackson_H

#include "ppModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class JohnsonJackson Declaration
\*---------------------------------------------------------------------------*/

class JohnsonJackson
:
    public ppModel
{

public:

    //- Runtime type information
    TypeName("JohnsonJackson");


    // Constructors

        //- Construct from components
        JohnsonJackson
        (
            const dictionary& ppDict,
            const phaseModel& phases,
            const phaseModel& phasew
        );


    //- Destructor
    virtual ~JohnsonJackson();

    // Member Functions

    tmp<volScalarField> pff
    (
        const volScalarField& alphas_,
        const volScalarField& alphasMinFriction,
        const dimensionedScalar& alphasMax,
        const dimensionedScalar& Fr,
        const dimensionedScalar& eta0,
        const dimensionedScalar& eta1
      ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
