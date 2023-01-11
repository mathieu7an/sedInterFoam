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

#include "partDynamicLagrangian.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void partDynamicLagrangian<BasicTurbulenceModel>::correctNut
(
    const tmp<volTensorField>& gradU
)
{
    this->nut_ = min((flm_/fmm_)*sqr(this->delta())*mag(dev(symm(gradU))),
                      nutMax_);
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}

template<class BasicTurbulenceModel>
void partDynamicLagrangian<BasicTurbulenceModel>::correctNut()
{
    correctNut(fvc::grad(this->U_));
}

template<class BasicTurbulenceModel>
volScalarField partDynamicLagrangian<BasicTurbulenceModel>::deltaStar
(
    const tmp<volScalarField>& delta,
    const tmp<volScalarField>& alpha,
    const tmp<volVectorField>& Ur,
    const tmp<volScalarField>& K
)
{
    const rhoField& rho = this->rho_;
    const volScalarField taup(rho/(max(K, Ksmall_)*(1-alpha)));
    const volScalarField deltaStarr(delta/(taup*mag(max(Ur, UrSmall_))));
    return deltaStarr;
}

template<class BasicTurbulenceModel>
volScalarField partDynamicLagrangian<BasicTurbulenceModel>::fDelta
(
    const tmp<volScalarField>& delta
)
{
    const volScalarField fDeltaa(sqr(delta)/(Cf_+sqr(delta)));
    return fDeltaa;
}

template<class BasicTurbulenceModel>
volScalarField partDynamicLagrangian<BasicTurbulenceModel>::hAlpha
(
    const tmp<volScalarField>& alpha
)
{
    return -tanh(alpha/Ch1_)*sqrt(alpha/alphasMax_)*sqr(1-alpha/alphasMax_)
            *(1-Ch2_*alpha/alphasMax_+Ch3_*sqr(alpha/alphasMax_));
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
partDynamicLagrangian<BasicTurbulenceModel>::partDynamicLagrangian
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    LESeddyViscosity<BasicTurbulenceModel>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),
    flm_
    (
        IOobject
        (
            IOobject::groupName("flm", this->U_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    fmm_
    (
        IOobject
        (
            IOobject::groupName("fmm", this->U_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    Cy_
    (
        IOobject
        (
            "Cy",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0)
    ),
    theta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "theta",
            this->coeffDict_,
            1.5
        )
    ),
    Cf_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cf",
            this->coeffDict_,
            dimensionSet(2, -6, 0, 0, 0),
            0.15
        )
    ),
    Ch1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ch1",
            this->coeffDict_,
            0.1
        )
    ),
    Ch2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ch2",
            this->coeffDict_,
            1.88
        )
    ),
    Ch3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ch3",
            this->coeffDict_,
            5.16
        )
    ),
    D_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "D",
            this->coeffDict_,
            0.0
        )
    ),
    alphasMax_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphasMax",
            this->coeffDict_,
            0.6
        )
    ),
    nutMax_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "nutaMax",
            this->coeffDict_,
            dimLength*dimVelocity,
            1e-3
        )
    ),
    Ksmall_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ksmall",
            this->coeffDict_,
            dimensionSet(1, -3, -1, 0, 0),
            VSMALL
        )
    ),
    UrSmall_
    (
        dimensioned<vector>::lookupOrAddToDict
        (
            "UrSmall",
            this->coeffDict_,
            dimensionSet(0, 1, -1, 0, 0),
            vector(SMALL, SMALL, SMALL)
        )
    ),
    CyMax_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CyMax",
            this->coeffDict_,
            1.0
        )
    ),
    onx_
    (
        IOobject
        (
            "onx",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        tensor(1, 0, 0, 0, 0, 0, 0, 0, 0)
    ),
    ony_
    (
        IOobject
        (
            "ony",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        tensor(0, 0, 0, 0, 1, 0, 0, 0, 0)
    ),
    onz_
    (
        IOobject
        (
            "onz",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        tensor(0, 0, 0, 0, 0, 0, 0, 0, 1)
    ),
    simpleFilter_(U.mesh()),
    filterPtr_(LESfilter::New(U.mesh(), this->coeffDict())),
    filter_(filterPtr_()),

    flm0_("flm0", flm_.dimensions(), 0.0),
    fmm0_("fmm0", fmm_.dimensions(), VSMALL)
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool partDynamicLagrangian<BasicTurbulenceModel>::read()
{
    if (LESeddyViscosity<BasicTurbulenceModel>::read())
    {
        filter_.read(this->coeffDict());
        theta_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
void partDynamicLagrangian<BasicTurbulenceModel>::correct()
{
    if (not this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    LESeddyViscosity<BasicTurbulenceModel>::correct();

    tmp<volTensorField> tgradU(fvc::grad(U));
    const volTensorField& gradU = tgradU();

    volSymmTensorField S(dev(symm(gradU)));
    volScalarField magS(mag(S));
    volScalarField magSqrS(magSqr(S));

    volVectorField Uf(filter_(U));
    volSymmTensorField Sf(dev(symm(fvc::grad(Uf))));
    volScalarField magSf(mag(Sf));
    volScalarField magSqrSf(magSqr(Sf));
    
    volSymmTensorField L(dev(filter_(sqr(U)) - (sqr(filter_(U)))));
    volSymmTensorField M
    (
        2.0*sqr(this->delta())*(filter_(magS*S) - 4.0*magSf*Sf)
    );

    volScalarField invT
    (
        alpha*rho*(1.0/(theta_.value()*this->delta()))
        *pow(flm_*fmm_, 1.0/8.0)
    );

    volScalarField LM(L && M);

    fvScalarMatrix flmEqn
    (
        fvm::ddt(alpha, rho, flm_)
      + fvm::div(alphaRhoPhi, flm_)
     ==
        invT*LM
      - fvm::Sp(invT, flm_)
      + fvOptions(alpha, rho, flm_)
    );
    flmEqn.relax();
    fvOptions.constrain(flmEqn);
    flmEqn.solve();
    fvOptions.correct(flm_);
    bound(flm_, flm0_);
    volScalarField MM(M && M);

    fvScalarMatrix fmmEqn
    (
        fvm::ddt(alpha, rho, fmm_)
      + fvm::div(alphaRhoPhi, fmm_)
     ==
        invT*MM
      - fvm::Sp(invT, fmm_)
      + fvOptions(alpha, rho, fmm_)
    );

    fmmEqn.relax();
    fvOptions.constrain(fmmEqn);
    fmmEqn.solve();
    fvOptions.correct(fmm_);
    bound(fmm_, fmm0_);

    correctNut(gradU);
    
    volScalarField Lkk(1.0/3.0*tr((filter_(sqr(U)))-sqr(filter_(U))));

    volScalarField Mkk
    (
        -2.0*sqr(this->delta())*(filter_(magSqrS) - 4.0*magSqrSf)
    );
    
    Lkk.correctBoundaryConditions();
    Mkk.correctBoundaryConditions();

    volScalarField Lkka(fvc::average(Lkk*Mkk));
    volScalarField Mkka(fvc::average(Mkk*Mkk));

    Cy_ = min(Lkka/Mkka, CyMax_);
    bound(Cy_, VSMALL);
}

template<class BasicTurbulenceModel>
tmp<volScalarField> partDynamicLagrangian<BasicTurbulenceModel>::
spherSigmaSGS()
{
    const volVectorField& U = this->U_;
    tmp<volTensorField> tgradU(fvc::grad(U));
    const volTensorField& gradU = tgradU();

    volSymmTensorField S(dev(symm(gradU)));
    volScalarField magSqrS(magSqr(S));

    const tmp<volScalarField> spherSigma(Cy_*sqr(this->delta())*magSqrS);

    return spherSigma;
}

template<class BasicTurbulenceModel>
tmp<volTensorField> partDynamicLagrangian<BasicTurbulenceModel>::
Ksgs()
{
    const alphaField& alpha = this->alpha_;
    const volVectorField& Us = this->U_;
    const volVectorField& Uw = this->mesh().objectRegistry::template
                               lookupObjectRef<volVectorField>("Uw");
    const volVectorField Ur(Uw - Us);
    const volScalarField& K = this->mesh().objectRegistry::template
                               lookupObjectRef<volScalarField>("K");
    volVectorField Usf(filter_(Us));
    volVectorField Uwf(filter_(Uw));
    volScalarField alphaf(filter_(alpha));
    const volVectorField Urf(Uwf - Usf);

    volVectorField gradAlpha(fvc::grad(alpha));
    tmp<volTensorField> tgradUs(fvc::grad(Us));
    const volTensorField& gradUs = tgradUs();

    volVectorField gradAlphaf(fvc::grad(alphaf));
    tmp<volTensorField> tgradUsf(fvc::grad(Usf));
    const volTensorField& gradUsf = tgradUsf();

    volScalarField deltaStarr(deltaStar(this->delta(), alpha, Ur, K));
    volScalarField deltaStarrf(deltaStar(2*this->delta(), alphaf, Urf, K));

    volVectorField L(filter_(alpha*Ur) - alphaf*Urf);
    volVectorField M
    (
        filter_(fDelta(deltaStarr)*hAlpha(alpha)*alpha*Ur)
      - fDelta(deltaStarrf)*hAlpha(alphaf)*alphaf*Urf
    );
    volVectorField H
    (
        4.0*D_*sqr(this->delta())*(gradAlphaf&gradUsf)
      - 1.0*D_*sqr(this->delta())*filter_(gradAlpha&gradUs)
    );

    L.correctBoundaryConditions();
    M.correctBoundaryConditions();

    volScalarField LHMx
    (
        fvc::average
        (
            (L.component(vector::X))//-H.component(vector::X))
            *M.component(vector::X)
        )
    );
    volScalarField LHMy
    (
        fvc::average
        (
            (L.component(vector::Y))//-H.component(vector::Y))
            *M.component(vector::Y)
        )
    );
    volScalarField LHMz
    (
        fvc::average
        (
            (L.component(vector::Z))//-H.component(vector::Z))
            *M.component(vector::Z)
        )
    );
    volScalarField MMx
    (
        fvc::average
        (
            M.component(vector::X)*M.component(vector::X)
        )
    );
    volScalarField MMy
    (
        fvc::average
        (
            M.component(vector::Y)*M.component(vector::Y)
        )
    );
    volScalarField MMz
    (
        fvc::average
        (
            M.component(vector::Z)*M.component(vector::Z)
        )
    );

    MMx.max(VSMALL);
    MMy.max(VSMALL);
    MMz.max(VSMALL);

    volScalarField Kxx(LHMx/MMx);
    volScalarField Kyy(LHMy/MMy);
    volScalarField Kzz(LHMz/MMz);

    const tmp<volTensorField> Ksgss(Kxx*onx_ + Kyy*ony_ + Kzz*onz_);

    return fDelta(deltaStarr)*hAlpha(alpha)*Ksgss;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
