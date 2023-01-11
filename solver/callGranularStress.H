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

gradUsT = fvc::grad(Us)().T();
gradUfT = fvc::grad(Uf)().T();
volScalarField magDtensor
(
    ::sqrt(2.0)*mag(symm(gradUsT))
);

if (granularRheology.dilatancy())
{
    //// Solving the evolution of plastic volume fraction

    fvScalarMatrix phi_pl_Eqn
    (
     fvm::ddt(alphasPlastic)
     + fvm::div(phis, alphasPlastic, "div(phis,alphasPlastic)")
     - fvm::Sp(fvc::div(phis), alphasPlastic)
     ==
     -fvc::Sp(delta*magDtensor, alphasPlastic) // explicit term
    );
    phi_pl_Eqn.relax();
    phi_pl_Eqn.solve();

    alphasPlastic.max(0.48);
    alphasPlastic.min(0.585);
}
else
{
    alphasPlastic=alphasMinFriction*Unity;
}

// Pff corresponds to "contact" pressure which prevent the granular phase from
// exceeding the maximum volume fraction.
pff = pp->pff(alphas, alphasPlastic, alphasMax, Fr_, eta0_, eta1_);

//
//
//  Kinetic THEORY PART
//
//
if (kineticTheory.on())
{
    dimensionedScalar Dsmall2
    (
        "Dsmall2",
        dimensionSet(0, 0, -2, 0, 0, 0, 0),
        1e-8
    );
    Dsmall2 = sqr(Dsmall);

    // Compute Kinetic Theory including granular temperature solution
    kineticTheory.solve
    (
        gradUsT, turbulencef->k(), turbulencef->epsilon(),
        turbulencef->nut(), B, runTime
    );

    // Compute Frictional viscosity
    volScalarField muEff_f
    (
        pff*Foam::sin(kineticTheory.phi())
       /sqrt(pow(magDtensor, 2) + Dsmall2)
    );

    // the actual expression for nuFrs depends on the way this term is
    // discretized in UsEqn
    nuFrs = muEff_f/rhos;

    // Compute nuEfff the total fluid phase viscosity
    nuEfff = nuf + turbulencef->nut();

    // Update solid phase viscosities and collisional pressure
    nuEffs = kineticTheory.mus()/((alphas + alphasSmall)*rhos)
             + turbulences->nut();
    lambdaUs = kineticTheory.lambda();
    ps = kineticTheory.ps();

    if (debugInfo)
    {
        Info<< "Contact pressure  pff: Min =" << gMin(pff)
            <<", Max =" << gMax(pff)<<endl;
        Info<< "Collisional press. ps: Min =" << gMin(ps)
            <<", Max =" << gMax(ps)<<endl;
    }
}

//
//
//  GRANULAR RHEOLOGY PART
//
//
else if (granularRheology.on())
{
// Solve granular rheology
    granularRheology.solve(magDtensor, pff, alphasSmall, Dsmall);

// Particulate pressure and viscosity
    ps = granularRheology.ps();

    delta=granularRheology.delta();

    if (debugInfo)
    {
        Info<< "Contact pressure  pff: Min =" << gMin(pff)
            <<", Max =" << gMax(pff)<<endl;
        Info<< "Shear ind. press.  ps: Min =" << gMin(ps)
            <<", Max =" << gMax(ps)<<endl;
    }
// the actual expression for nuFrs depends on the way this term is discretized
// in UsEqn
    nuFrs = granularRheology.mus()/ rhos;

    nuEffs = sqr(Ct)*turbulencef->nut() + nus + turbulences->nut();
    lambdaUs = nus*rhos*scalar(0.0);

// Compute nuEfff the total fluid phase viscosity
    nuEfff = turbulencef->nut() + granularRheology.nuvf();
}
//
//
//  If not using kinetic theory or granular rheology Ct model is used
//
//
else
{
    nuEffs = sqr(Ct)*turbulencef->nut() + nus;
    nuEfff = turbulencef->nut() + nuf;

    // set ps and nuFrs to zero if kineticTheory and granularRheology are off
    ps = ps*scalar(0.0);
    nuFrs = nus*scalar(0.0);
    nuEffs = nus + turbulences->nut();
    lambdaUs = nus*rhos*scalar(0.0);
}
// Add a numerical viscosity to damp the instabilities close to the outlet
if (spongeLayer)
{
    volScalarField XX(mesh.C().component(vector::X));
    nuEffs +=
    (
        pos(XX-xSmin)*nus*1e3
       *Foam::exp(-(xSmax-XX)/max(XX-xSmin, 1e-10*(xSmax-xSmin)))
    );
    nuEfff +=
    (
        pos(XX-xSmin)*nuf*1e3
       *Foam::exp(-(xSmax-XX)/max(XX-xSmin, 1e-10*(xSmax-xSmin)))
);
}
// Limit viscosities for numerical stability
nuFrs.min(nuMax);
nuEffs.min(nuMax);
nuEfff.min(nuMax);
