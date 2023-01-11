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

#include "granularRheologyModel.H"
#include "FrictionModel.H"
#include "PPressureModel.H"
#include "FluidViscosityModel.H"
#include "surfaceInterpolate.H"
#include "mathematicalConstants.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::granularRheologyModel::granularRheologyModel
(
    const Foam::phaseModel& phases,
    const Foam::volScalarField& rhof,
    const Foam::volScalarField& nuf,
    const Foam::volScalarField& ps,
    const Foam::dimensionedScalar& Dsmall
)
:
    alphas_(phases.alpha()),
    phis_(phases.phi()),
    rhos_(phases.rho()),
    ds_(phases.d()),
    rhof_(rhof),
    nuf_(nuf),

    ps_new_value(ps),

    granularRheologyProperties_
    (
        IOobject
        (
            "granularRheologyProperties",
            alphas_.time().constant(),
            alphas_.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    granularRheology_
    (
        granularRheologyProperties_.get<Switch>("granularRheology")
    ),
    granularDilatancy_
    (
        granularRheologyProperties_.get<Switch>("granularDilatancy")
    ),
    FrictionModel_
    (
        granularRheologyModels::FrictionModel::New
        (
            granularRheologyProperties_
        )
    ),
    PPressureModel_
    (
        granularRheologyModels::PPressureModel::New
        (
            granularRheologyProperties_
        )
    ),
    FluidViscosityModel_
    (
        granularRheologyModels::FluidViscosityModel::New
        (
            granularRheologyProperties_
        )
    ),
    alphasMaxG_
    (
        granularRheologyProperties_.getOrDefault
        (
            "alphasMaxG",
            dimensionedScalar("alphasMaxG",
                          dimensionSet(0, 0, 0, 0, 0, 0, 0),
                          0.6)
        )
    ),
    muss_
    (
        granularRheologyProperties_.getOrDefault
        (
            "muss",
            dimensionedScalar("muss",
                          dimensionSet(0, 0, 0, 0, 0, 0, 0),
                          0.38)
        )
    ),
    mu2_
    (
        granularRheologyProperties_.getOrDefault
        (
            "mu2",
            dimensionedScalar("mu2",
                          dimensionSet(0, 0, 0, 0, 0, 0, 0),
                          0.64)
        )
    ),
    I0_
    (
        granularRheologyProperties_.getOrDefault
        (
            "I0",
            dimensionedScalar("I0",
                          dimensionSet(0, 0, 0, 0, 0, 0, 0),
                          0.3)
        )
    ),
    Bphi_
    (
        granularRheologyProperties_.getOrDefault
        (
            "Bphi",
            dimensionedScalar("Bphi",
                          dimensionSet(0, 0, 0, 0, 0, 0, 0),
                          0.31)
        )
    ),
    n_
    (
        granularRheologyProperties_.getOrDefault
        (
            "n",
            dimensionedScalar("n",
                          dimensionSet(0, 0, 0, 0, 0, 0, 0),
                          2.5)
        )
    ),
    BulkFactor_
    (
        granularRheologyProperties_.getOrDefault
        (
                "BulkFactor",
                dimensionedScalar("BulkFactor",
                    dimensionSet(0, 0, 0, 0, 0, 0, 0),
                    0)
        )
    ),
    alphas_c_
    (
        granularRheologyProperties_.getOrDefault
        (
                "alphas_c",
                dimensionedScalar("alphas_c",
                    dimensionSet(0, 0, 0, 0, 0, 0, 0),
                    0.585)
        )
    ),

    K_dila_
    (
        granularRheologyProperties_.getOrDefault
        (
                "K_dila",
                dimensionedScalar("K_dila",
                    dimensionSet(0, 0, 0, 0, 0, 0, 0),
                    0)
        )
    ),
    relaxPs_
    (
        granularRheologyProperties_.getOrDefault
        (
            "relaxPs",
            dimensionedScalar("relaxPs",
                          dimensionSet(0, 0, 0, 0, 0, 0, 0),
                          1)
        )
    ),
    PsMin_
    (
        granularRheologyProperties_.getOrDefault
        (
            "PsMin",
            dimensionedScalar("PsMin",
                          dimensionSet(0, 1, 0, 0, 0, 0, 0),
                          1e-6)
        )
    ),
    tau_inv_min_
    (
        granularRheologyProperties_.getOrDefault
        (
            "tau_inv_min",
            dimensionedScalar("tau_inv_min",
                          dimensionSet(0, 0, -1, 0, 0, 0, 0),
                          1e-12)
        )
    ),
    muI_
    (
        IOobject
        (
            "muI",
            alphas_.time().timeName(),
            alphas_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alphas_.mesh(),
        dimensionedScalar("zero", alphas_.dimensions(), 0.0)
    ),
    mus_
    (
        IOobject
        (
            "mus",
            alphas_.time().timeName(),
            alphas_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alphas_.mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -1, 0, 0), 0.0)
    ),
    lambda_
    (
        IOobject
        (
            "lambda",
            alphas_.time().timeName(),
            alphas_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alphas_.mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -1, 0, 0), 0.0)
    ),

    ps_
    (
        IOobject
        (
            "ps",
            alphas_.time().timeName(),
            alphas_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alphas_.mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -2, 0, 0), 0.0)
    ),
    p_p_total_
    (
        IOobject
        (
            "p_p_total",
            alphas_.time().timeName(),
            alphas_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alphas_.mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -2, 0, 0), 0.0)
    ),

    delta_
    (
        IOobject
        (
            "delta",
            alphas_.time().timeName(),
            alphas_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alphas_.mesh(),
        dimensionedScalar("zero", alphas_.dimensions(), 0.0)
    ),

    nuvf_
    (
        IOobject
        (
            "nuvf",
            alphas_.time().timeName(),
            alphas_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alphas_.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 2, -1, 0, 0), 0.0)
    ),
    I_
    (
        IOobject
        (
            "I",
            alphas_.time().timeName(),
            alphas_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alphas_.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::granularRheologyModel::~granularRheologyModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::granularRheologyModel::solve
(
    const volScalarField& magD,
    const volScalarField& pf,
    const dimensionedScalar& alphasSmall,
    const dimensionedScalar& Dsmall
)
{
    if (granularRheology_ == false)
    {
        return;
    }
    dimensionedScalar Dsmall2
    (
        "Dsmall2",
        dimensionSet(0, 0, -2, 0, 0, 0, 0),
        1e-8
    );
    Dsmall2 = sqr(Dsmall);
    //
    // compute the particulate velocity shear rate
    //
    volScalarField magD2(pow(magD, 2));

    //
    // Shear induced particulate pressure
    //
    ps_ = PPressureModel_->ps
    (
        pf, Bphi_, rhos_, ds_, rhof_, nuf_, magD,
        alphas_, alphasMaxG_, alphasSmall
    );

    // Relaxing shear induced particulate pressure
    //  relaxPs_ controls the relaxation of ps. Low values lead to relaxed ps
    //  whereas large value are prone to numerical error
    volScalarField tau_inv_par(relaxPs_*alphas_*magD);
    tau_inv_par.max(tau_inv_min_);

    fvScalarMatrix psEqn
    (
         fvm::ddt(ps_new_value)
         + fvm::div(phis_, ps_new_value, "div(phis,ps_new_value)")
         - fvm::Sp(fvc::div(phis_), ps_new_value)
        ==
        tau_inv_par*(ps_)
        -fvm::Sp(tau_inv_par, ps_new_value)
    );
    psEqn.relax();
    psEqn.solve();

    ps_new_value.max(0.0);

    ps_=ps_new_value;
    ps_.correctBoundaryConditions();
    //total particle pressure(shear induced+contact contributions)
    p_p_total_ = mag(ps_new_value+pf);
    p_p_total_.max(PsMin_);

    //  Compute the particulate friction coefficient
    muI_ = FrictionModel_->muI(muss_, mu2_, I0_, p_p_total_, rhos_, ds_, rhof_,
                               nuf_, magD, Dsmall);

    //  Compute the inertial/viscous number
    I_ = FrictionModel_->I(p_p_total_, rhos_, ds_, rhof_, nuf_, magD);

    // Dilatancy model
    if (granularDilatancy_)
    {
    //delta_ = DilatancyModel_->delta(K_dila_, alphas_c_, alphas_, magD,
    // ds_, rhof_, nuf_, p_p_total_, PsMin);
        volScalarField alphasEq_
        (
            PPressureModel_->alphasEq
            (
                p_p_total_,
                Bphi_,
                rhos_,
                ds_,
                rhof_,
                nuf_,
                magD,
                alphas_c_
            )
        );
        delta_ = K_dila_*(alphas_ - alphasEq_);

        delta_.min( 0.5);
        delta_.max(-0.5);
    }
    //  Compute the regularized particulate viscosity
    mus_ = muI_* p_p_total_ / pow(magD2 + Dsmall2, 0.5);

    // Compute bulk viscosity (by default BulkFactor = 0)s
    lambda_ = BulkFactor_*p_p_total_ / pow(magD2 + Dsmall2, 0.5);

    // Compute the Effective fluid viscosity
    nuvf_ = FluidViscosityModel_->nuvf(alphas_, nuf_, alphasMaxG_, alphasSmall,
                                       n_);
}
// ************************************************************************* //
