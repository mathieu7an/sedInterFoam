/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
scalar pi(constant::mathematical::pi);

surfaceScalarField rhoff(fvc::interpolate(rhof));

// following is general form of mean flow + sine wave and 2nd order stokes waves
// Add oscilitary pressure gradient
gradPOSC =
(
    gradPMEAN
  + gradPAMP1*Foam::cos(initTheta1*pi/180.0 + 1.0*2.0*pi*runTime.value()/oscpT)
  + gradPAMP2*Foam::cos(initTheta2*pi/180.0 + 2.0*2.0*pi*runTime.value()/oscpT)
  + gradPAMP3*Foam::cos(initTheta3*pi/180.0 + 3.0*2.0*pi*runTime.value()/oscpT)
  + gradPAMP4*Foam::cos(initTheta4*pi/180.0 + 4.0*2.0*pi*runTime.value()/oscpT)
  + gradPAMP5*Foam::cos(initTheta5*pi/180.0 + 5.0*2.0*pi*runTime.value()/oscpT)
);
// Compute mixture density
rho = alphas*rhos+alphaf*rhof;
if (SUSaniso)
{
    suss = SUS*(K/rhos)*alphaf/(alphas+alphasSmall)*turbulencef->nut()
                   *((SUS_I*iso-SUS_A*aniso) & fvc::grad(alphas));
    susf = SUS*(K/rhof)*turbulencef->nut()
                   *((SUS_I*iso-SUS_A*aniso) & fvc::grad(alphas));
}
// limiter for the explicit stress terms to avoid numerical instability and
// too strong limitation on the time step
volScalarField nuFrsExp(nuFrs);
nuFrsExp.min(nuMaxExp);
nuFrsExp.correctBoundaryConditions();
// Compute explicit stress term for phase s
volTensorField Rcs
(
    "Rcs",
    ((2.0/3.0)*I)*((nuEffs+nuFrsExp/(alphas+alphasSmall))*tr(gradUsT)
    + spherSigmaSGSs)
    - (nuEffs+nuFrsExp/(alphas+alphasSmall))*gradUsT
);
// Compute concentration gradient stress contribution on cell faces for phase s
surfaceScalarField phiRs
(
   -fvc::interpolate(nuEffs)*mesh.magSf()*(fvc::snGrad(alphas))
   /fvc::interpolate(alphas+ alphasSmall)
);
// Assemble momentum balance equation for phase s
fvVectorMatrix UsEqn
(
    (scalar(1.0) + Cvm*rhof*alphaf/rhos)*
    (
        fvm::ddt(rhos, Us)
      + fvm::div(rhos*phis, Us)
      - fvm::Sp(rhos*fvc::div(phis), Us)
    )
// effective stress term from dispersed phase
  - scalar(1.0)/(alphas+alphasSmall)*fvm::laplacian(rhos*nuFrs, Us)
  - fvm::laplacian(rhos*nuEffs, Us)
  + fvm::div(rhos*phiRs, Us, "div(phiRs,Us)")
  - fvm::Sp(rhos*fvc::div(phiRs), Us)
  + fvc::div(rhos*Rcs)
  + rhos*(fvc::grad(alphas)/(alphas+alphasSmall) & Rcs)
    ==
//  - gradPf                  // contact pressure gradient moved to p-equation
  - rhos*suss                      // anistropic turbulent suspension term
  - fvc::grad(pff)/(alphas + alphasSmall)
  - alphaf*(liftCoeff - Cvm*rhof*DDtUf)  // lift and added mass forces
  //+ fvOptions(rhos*alphas, Us)
);
UsEqn.relax();
UsEqn += fvm::Sp(alphaf*K, Us);
fvOptions.constrain(UsEqn);
Us.correctBoundaryConditions();
//fvOptions.correct(Us);
// Compute explicit stress term for phase b
volTensorField Rcf
(
    "Rcf",
    ((2.0/3.0)*I)*(nuEfff*tr(gradUfT) + spherSigmaSGSf)
    - nuEfff*gradUfT
);
// Compute concentration gradient stress contribution on cell faces for phasew
surfaceScalarField phiRf
(
  - fvc::interpolate(rhof*nuEfff)*mesh.magSf()*(fvc::snGrad(alphaf))
   /fvc::interpolate(alphaf)
);

// Assemble momentum balance equation for phase w
fvVectorMatrix UfEqn
(
    (scalar(1.0) + Cvm*rhof*alphas/rhof)*
    (
        fvm::ddt(rhof,Uf)
      + fvm::div(rhoPhi, Uf)
      - fvm::Sp(rhof*fvc::div(phif), Uf)
    )
  - fvm::laplacian(rhof*nuEfff, Uf) 
  + fvm::div(phiRf, Uf, "div(phiRf,Uf)")
  - fvm::Sp(fvc::div(phiRf), Uf)
  + fvc::div(rhof*Rcf)
  + (rhof*fvc::grad(alphaf)/(alphaf) & Rcf)
    ==
    rhof*susf     // anistropic turbulent suspension term
  //+ sigma*Kg*fvc::grad(gamma)
  //- gh*fvc::grad(rhof)
  + alphas*(liftCoeff + Cvm*rhof*DDtUs) // lift & added mass forces
  //+ fvOptions(rhof*alphaf, Uf)
);
UfEqn.relax();
UfEqn += fvm::Sp(alphas*K, Uf);
fvOptions.constrain(UfEqn);
Uf.correctBoundaryConditions();
//fvOptions.correct(Uf);