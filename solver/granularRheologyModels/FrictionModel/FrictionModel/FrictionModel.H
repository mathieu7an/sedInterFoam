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

Class
    Foam::FrictionModel

Description

    This file is in the \dir FrictionModel/FrictionModel subdirectory
    inside the \dir FrictionModel subdirectory

SourceFiles
    FrictionModel.C

\*---------------------------------------------------------------------------*/

#ifndef FrictionModel_H
#define FrictionModel_H

#include "dictionary.H"
#include "volFields.H"
#include "dimensionedTypes.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace granularRheologyModels
{

/*---------------------------------------------------------------------------*\
                           Class FrictionModel Declaration
\*---------------------------------------------------------------------------*/

class FrictionModel
{
    // Private member functions

    //- Disallow default bitwise copy construct
    FrictionModel(const FrictionModel&);

    //- Disallow default bitwise assignment
    void operator=(const FrictionModel&);


protected:

    // Protected data

    const dictionary& dict_;


public:

    //- Runtime type information
    TypeName("FrictionModel");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        FrictionModel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );


    // Constructors

    //- Construct from components
    FrictionModel(const dictionary& dict);


    // Selectors

    static autoPtr<FrictionModel> New
    (
        const dictionary& dict
    );


    //- Destructor
    virtual ~FrictionModel();


    // Member Functions

    virtual tmp<volScalarField> muI
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
    ) const = 0;

    virtual tmp<volScalarField> I
    (
        const volScalarField& ps,
        const dimensionedScalar& rhos,
        const dimensionedScalar& ds,
        const volScalarField& rhof,
        const volScalarField& nuf,
        const volScalarField& magD
    ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace granularRheologyModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
