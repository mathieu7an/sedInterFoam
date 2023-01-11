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
{
    word schemeA("div(alphas,nu)");

    // Compute volume fraction of phase s and f at cell faces
    surfaceScalarField alphasf(fvc::interpolate(alphas));
    surfaceScalarField alphaff(scalar(1.0) - alphasf);
    surfaceScalarField rhoff(fvc::interpolate(rhof));

    // Compute inverse matrix for momentum balance of phase s and f
    //     at cell centers
    volScalarField rUsA(1.0/UsEqn.A());
    volScalarField rUfA(1.0/UfEqn.A());

    // Interpolate inverse matrix for momentum balance of phase s and f
    //     at cell faces
    surfaceScalarField rUsAf(fvc::interpolate(rUsA));
    surfaceScalarField rUfAf(fvc::interpolate(rUfA));

    // limiter for the explicit stress terms to avoid numerical instability and
    // too strong limitation on the time step
    volScalarField nuFrsExp(nuFrs);
    nuFrsExp.min(nuMaxExp);
    nuFrsExp.correctBoundaryConditions();

    // interpolate 1/(alpha+alphasSmall) on the faces
    surfaceScalarField oneOverAlphas
    (
        fvc::interpolate
        (
            scalar(1.0)/(alphas+alphasSmall),
            "interpolate(oneoveralphas)"
        )
    );

    // local Schmidt number
    if (SUSlocal)
    {
        volScalarField ustar
        (
            sqrt(mag(alphaf*nuEfff*(fvc::grad(Uf)+gradUfT)))
        );
        SUStilde = scalar(1.0) + scalar(2.0)*pow(Ufall/(ustar+Usmall), 2);
        SUStilde.min(scalar(3.0));
        SUS = SUStilde;

        Info<<"ustarmax"<<max(ustar).value()<<endl;
    }

    // PIMPLE algorithm loop
    while (pimple.correct())
    {
        // Compute predicted velocity of phase a and b at cell centers
        volVectorField HsbyA(Us);
        HsbyA = rUsA*UsEqn.H();
        volVectorField HfbyA(Uf);
        HfbyA = rUfA*UfEqn.H();

        // Momentum flux at the cell faces for phase a
        surfaceScalarField phiDrags
        (
          // grad(pa)
          - (oneOverAlphas/rhos*fvc::snGrad(ps))*rUsAf*mesh.magSf()
          // grad(Pf)
          - (oneOverAlphas/rhos*fvc::snGrad(pff))*rUsAf*mesh.magSf()
          // diagonal part of stress
          + (oneOverAlphas/rhos*fvc::snGrad(lambdaUs*tr(gradUsT)))
           *rUsAf*mesh.magSf()
          // diagonal part of stress 2
          - (
                scalar(2.0/3.0)*oneOverAlphas
               *fvc::snGrad((alphas+alphasSmall)*(sqr(Ct)*turbulencef->k()
              + nuEffs*tr(gradUsT) + spherSigmaSGSs))
            )
           *rUsAf*mesh.magSf()
          - (scalar(2.0/3.0)*oneOverAlphas*fvc::snGrad(nuFrsExp*tr(gradUsT)))
           *rUsAf*mesh.magSf()
          // stress part 3
          + rUsAf*(oneOverAlphas
           *fvc::interpolate
            (
                fvc::div((alphas+alphasSmall)*nuEffs*gradUsT, schemeA)
            )
          & mesh.Sf())
          + rUsAf*(oneOverAlphas
           *fvc::interpolate(fvc::div(nuFrsExp*gradUsT, schemeA)) & mesh.Sf())
          // gravity
          + rUsAf
           *(
               (1.0-rhoff/rhos)*(gravity & mesh.Sf())
            )
          // mean pressure gradient
          + rUsAf*((alphaff*gradPOSC/(alphaff*(tilt*rhoff+(1.0-tilt)*rhos))) & mesh.Sf())
        );
    if (not SUSaniso)
    {
        // turbulent suspension term
        phiDrags += -(fvc::interpolate(SUS*(K/rhos)*turbulencef->nut()
                   *(alphaf/(alphas + alphasSmall)))
                   *fvc::snGrad(alphas))*rUsAf*mesh.magSf();
    }

        // Momentum flux at the cell faces for phase f
        surfaceScalarField phiDragf
        (
        // drag modification to account for particle fluctuations
        // -phiTempb
            // diagnoal part of stress
          - (fvc::interpolate(scalar(2.0)/(alphaf*scalar(3.0)))
           *fvc::snGrad(((nuEfff*tr(gradUfT) + spherSigmaSGSf))*alphaf))
           *rUfAf*mesh.magSf()
          // stress part 3
          + rUfAf*(fvc::interpolate(scalar(1.0)/alphaf)
           *fvc::interpolate
            (
                fvc::div(alphaf*nuEfff*gradUfT, schemeA)) & mesh.Sf()
            )
          // driving force
          + rUfAf*((gradPOSC*alphaff/(rhoff*alphaff)) & mesh.Sf())
          // Convective acceleration for 1D non-uniform simulation
          - rUfAf*fvc::interpolate(neg(alphas-0.3))*((alphaff*gradConvAcc/(alphaff*rhoff))
            & mesh.Sf())
        );
    if (not SUSaniso)
    {
        // turbulent suspension term
        phiDragf += fvc::interpolate(SUS*(K/rhof)
                        *turbulencef->nut())*fvc::snGrad(alphas)
                        *rUfAf*mesh.magSf();
    }
        //    fvOptions.makeRelative(phiHabyA);
        //    fvOptions.makeRelative(phiHbbyA);


        // Fix for gravity on outlet boundary.
        forAll(p_rbgh.boundaryField(), patchi)
        {
             if (isA<zeroGradientFvPatchScalarField>
                 (p_rbgh.boundaryField()[patchi]))
             {
                 phiDrags.boundaryFieldRef()[patchi] = 0.0;
                 phiDragf.boundaryFieldRef()[patchi] = 0.0;
             }
             if (isA<symmetryFvPatchScalarField>
                 (p_rbgh.boundaryField()[patchi]))
             {
                 phiDrags.boundaryFieldRef()[patchi] = 0.0;
                 phiDragf.boundaryFieldRef()[patchi] = 0.0;
             }
             if (p_rbgh.boundaryField().types()[patchi] ==
                 "symmetryPlane")
             {
                 phiDrags.boundaryFieldRef()[patchi] = 0.0;
                 phiDragf.boundaryFieldRef()[patchi] = 0.0;
             }
             if (p_rbgh.boundaryField().types()[patchi] ==
                 "groovyBC")
             {
                 phiDrags.boundaryFieldRef()[patchi] = 0.0;
                 phiDragf.boundaryFieldRef()[patchi] = 0.0;
             }
             if (p_rbgh.boundaryField().types()[patchi] ==
                 "outletInlet")
             {
                 phiDrags.boundaryFieldRef()[patchi] = 0.0;
                 phiDragf.boundaryFieldRef()[patchi] = 0.0;
             }
             if (p_rbgh.boundaryField().types()[patchi] ==
                 "fixedFluxPressure")
             {
                 phiDrags.boundaryFieldRef()[patchi] = 0.0;
                 phiDragf.boundaryFieldRef()[patchi] = 0.0;
             }
        }
        // Compute corrected velocity fluxes at cell faces
        //      (account for time evolution)
        surfaceScalarField phiUs
        (
            (fvc::interpolate(HsbyA) & mesh.Sf()) + rUsAf*fvc::ddtCorr(Us, phis)
        );

        surfaceScalarField phiUf
        (
            (fvc::interpolate(HfbyA) & mesh.Sf()) + rUfAf*fvc::ddtCorr(Uf, phif)
        );
        // Compute predicted velocity fluxes at cell faces
        //      (account for momentum fluxes computed at cell faces)
        surfaceScalarField phisHbyA(phiUs + phiDrags);
        surfaceScalarField phifHbyA(phiUf + phiDragf);

        // Face-drag coefficients
        surfaceScalarField rAUKds(fvc::interpolate(alphaf/rhos*K*rUsA));
        surfaceScalarField rAUKdf(fvc::interpolate(alphas/rhof*K*rUfA));
        // Compute predicted velocity fluxes at cell faces including
        //       including semi-implicit drag contribution
        surfaceScalarField phiHbyA
        (
            alphasf*(phisHbyA + rAUKds*phif)
          + alphaff*(phifHbyA + rAUKdf*phis)
        );

        // adjustPhi(phiHbyA, U, p_rbgh);

        // Assemble Poisson operator for the pressure
        surfaceScalarField Dp
        (
            "(rho*(1|A(U)))",
            alphasf*rUsAf/rhos + alphaff*rUfAf/rhoff
        );

        // Update the fixedFluxPressure BCs to ensure flux consistency
        setSnGrad<fixedFluxPressureFvPatchScalarField>
        (
            p_rbgh.boundaryFieldRef(),
            (
                phiHbyA.boundaryField()
              - alphasf.boundaryField()
               *(mesh.Sf().boundaryField() & Us.boundaryField())
              + alphaff.boundaryField()
               *(mesh.Sf().boundaryField() & Uf.boundaryField())
             )
            /(mesh.magSf().boundaryField()*Dp.boundaryField())
        );

        while (pimple.correctNonOrthogonal())
        {
            // Solve for the Poisson equation on the pressure
            fvScalarMatrix pEqn
            (
                fvm::laplacian(Dp, p_rbgh) == fvc::div(phiHbyA)
            );
            // Set pressure reference
            pEqn.setReference(pRefCell, getRefCellValue(p_rbgh, pRefCell));

#if (defined(OFVERSION) && (OFVERSION >= 700))
            pEqn.solve();
#else
            pEqn.solve(mesh.solver(p_rbgh.select(pimple.finalInnerIter())));
#endif

            if (pimple.finalNonOrthogonalIter())
            {
                // Compute pressure gradient at cell faces
                surfaceScalarField SfGradp(pEqn.flux()/Dp);
                // Compute corrected velocity fluxes at cell faces
                //       due to pressure gradient
                surfaceScalarField phiss(phisHbyA - rUsAf*SfGradp/rhos);
                surfaceScalarField phifs(phifHbyA - rUfAf*SfGradp/rhoff);
                // Compute mixture velocity flux at cell faces
                phi = phiHbyA - pEqn.flux();
                // Compute corrected relative velocity fluxes due to drag force
                //       at cell faces
                surfaceScalarField phir
                (
                    ((phiss + rAUKds*phifs) - (phifs + rAUKdf*phiss))
                   /(1 - rAUKds*rAUKdf)
                );
                // Compute corrected phase velocity fluxes at cell faces
                //       using updated semi-implicit drag
                phis = phi + alphaff*(phir);
                phif = phi - alphasf*(phir);

                // Compute contibuity error
                #include "continuityErrs.H"
                p_rbgh.relax();
                // Is it useful?
                SfGradp = pEqn.flux()/Dp;

                // Velocity correction step at cell centers
                volVectorField Ufs
                (
                    HfbyA + rUfA*fvc::reconstruct(phiDragf/rUfAf-SfGradp/rhoff)
                );
                volVectorField Uss
                (
                    HsbyA + rUsA*fvc::reconstruct(phiDrags/rUsAf-SfGradp/rhos)
                );
                // Compute corrected relative velocity due to drag force
                //       at cell centers
                volScalarField Ds(alphaf/rhos*K*rUsA);
                volScalarField Df(alphas/rhof*K*rUfA);
                volVectorField Ur(((1 - Df)*Uss - (1 - Ds)*Ufs)/(1 - Ds*Df));
                // Compute mixture velocity flux at cell centers
                U = alphas*(Uss + Ds*Uf) + alphaf*(Ufs + Df*Us);
                // Compute corrected phase velocity at cell centers
                //       using updated semi-implicit drag
                Us = U + alphaf*Ur;
                Us.correctBoundaryConditions();
                fvOptions.correct(Us);

                Uf = U - alphas*Ur;
                Uf.correctBoundaryConditions();
                fvOptions.correct(Uf);

                // Compute total pressure including hydrostatic contribution
                p = p_rbgh + rhof*gh;
                if (p_rbgh.needReference())
                {
                    p += dimensionedScalar
                    (
                        "p",
                        p.dimensions(),
                        pRefValue - getRefCellValue(p, pRefCell)
                    );
                    p_rbgh = p - rhof*gh;
                    p_rbgh.correctBoundaryConditions();
                }
                // Apply velocity limiters and correct boundary conditions
                #include "UsUfLimiter.H"
                Us.correctBoundaryConditions();
                // Update mixture velocity
                U = alphas*Us + alphaf*Uf;
                U.correctBoundaryConditions();
                Uw = Uf;
                Uw.correctBoundaryConditions();
                Ug = Uf;
                Ug.correctBoundaryConditions();
            }
        }
    }
}
