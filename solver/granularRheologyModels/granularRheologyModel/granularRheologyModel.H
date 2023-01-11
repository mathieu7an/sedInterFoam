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
    Foam::granularRheologyModel

Description
    Dense granular-flow rheology model for particle phase stress

    References:
    \verbatim
        GDRmidi (2004).
        On dense granular flows.
        The European Physical Journal E, 14, 341-365

        Forterre, Y. and Pouliquen, O. (2008).
        Flows of Dense Granular Media.
        Annual Review of Fluid Mechanics, 40, 1-24
    \endverbatim

    This file is in the \dir granularRheologyModels/granularRheologyModel
    subdirectory
    inside the \dir granularRheologyModels subdirectory

SourceFiles
    granularRheologyModel.C

\*---------------------------------------------------------------------------*/

#ifndef granularRheologyModel_H
#define granularRheologyModel_H

#include "dragModel.H"
#include "phaseModel.H"
#include "autoPtr.H"
#include "FrictionModel.H"
#include "PPressureModel.H"
#include "FluidViscosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class granularRheologyModel Declaration
\*---------------------------------------------------------------------------*/

class granularRheologyModel
{
    // Private data

    const volScalarField& alphas_;
    const surfaceScalarField& phis_;
    const dimensionedScalar& rhos_;
    const dimensionedScalar& ds_;
    const volScalarField& rhof_;
    const volScalarField& nuf_;

    volScalarField ps_new_value;

    //- dictionary holding the modeling info
    IOdictionary granularRheologyProperties_;

    //- use granular Rheology or not.
    Switch granularRheology_;

    //- use hydrostatic particulate pressure.
    Switch granularHydrostaticPressure_;

    //- use granular mixture model.
    Switch granularMixtureModel_;

    //- use dilatancy model or not.
    Switch granularDilatancy_;

    autoPtr<granularRheologyModels::FrictionModel> FrictionModel_;
    autoPtr<granularRheologyModels::PPressureModel> PPressureModel_;
    autoPtr<granularRheologyModels::FluidViscosityModel> FluidViscosityModel_;

    //- maximum packing
    const dimensionedScalar alphasMaxG_;

    //- Static friction coefficient
    const dimensionedScalar muss_;

    //- Dynamic friction coefficient
    const dimensionedScalar mu2_;

    //- Empirical rheology coefficient
    const dimensionedScalar I0_;

    //- Empirical Dilatancy coefficient
    const dimensionedScalar Bphi_;

    //- Exponent of the Krieger-Dougherty viscosity
    const dimensionedScalar n_;

    //- Bulk viscosity amplification factor
    const dimensionedScalar BulkFactor_;

    // critical volume fraction
    const dimensionedScalar alphas_c_;

    // dilatancy factor
    const dimensionedScalar K_dila_;

    //- Relaxation factor for pa
    const dimensionedScalar relaxPs_;

    //- Minimum particle pressure (regularisation)
    const dimensionedScalar PsMin_;

    //- Minimum relaxation timescale
    const dimensionedScalar tau_inv_min_;

    //- The friction coefficient
    volScalarField muI_;

    //- The granular viscosity
    volScalarField mus_;

    //- The granular bulk viscosity
    volScalarField lambda_;

    //- The granular pressure
    volScalarField ps_;

    //- Total particle  pressure
    volScalarField p_p_total_;

    ////- The  dilatancy angle
    volScalarField delta_;

    //- Fluid Effective viscosity (viscous contribution)
    volScalarField nuvf_;

    //-Viscous number
    volScalarField I_;

    // Private Member Functions

    //- Disallow default bitwise copy construct
    granularRheologyModel(const granularRheologyModel&);

    //- Disallow default bitwise assignment
    void operator=(const granularRheologyModel&);


public:

    // Constructors

    //- Construct from components
    granularRheologyModel
    (
        const phaseModel& phases,
        const volScalarField& rhof,
        const volScalarField& nuf,
        const volScalarField& p_p_total,
        const dimensionedScalar& Dsmall
    );


    //- Destructor
    virtual ~granularRheologyModel();


    // Member Functions

    void solve
    (
        const volScalarField& magD,
        const volScalarField& pf,
        const dimensionedScalar& alphasSmall,
        const dimensionedScalar& Dsmall
    );

    bool on() const
    {
        return granularRheology_;
    }

    bool mixture() const
    {
        return granularMixtureModel_;
    }
    bool dilatancy() const
    {
        return granularDilatancy_;
    }

    const volScalarField& mus() const
    {
        return mus_;
    }

    const volScalarField& ps() const
    {
        return ps_;
    }

    const dimensionedScalar& alphasMaxG() const
    {
        return alphasMaxG_;
    }
    const volScalarField& p_p_total() const
    {
        return p_p_total_;
    }
    const volScalarField& I() const
    {
        return I_;
    }
    const volScalarField& muI() const
    {
        return muI_;
    }
    const volScalarField& delta() const
    {
        return delta_;
    }
    const volScalarField& lambda() const
    {
        return lambda_;
    }
    const volScalarField& nuvf() const
    {
        return nuvf_;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
