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
    // Compute volume fraction of phase a and b at cell faces
    surfaceScalarField alphasf(fvc::interpolate(alphas));
    surfaceScalarField alphaff(scalar(1.0) - alphasf);

    // Compute inverse matrix for momentum balance of phase a and b
    //     at cell centers
    volScalarField rUsA(1.0/UsEqn.A());
    volScalarField rUfA(1.0/UfEqn.A());

    // Interpolate inverse matrix for momentum balance of phase a and b
    //     at cell faces
    surfaceScalarField rUsAf(fvc::interpolate(rUsA));
    surfaceScalarField rUfAf(fvc::interpolate(rUfA));

    // interpolate 1/(alpha+alphasSmall) on the faces
    surfaceScalarField oneOverAlphas
    (
        fvc::interpolate
        (
            scalar(1.0)/(alphas+alphasSmall),
            "interpolate(oneoveralphas)"
        )
    );

    surfaceScalarField oneOverAlphaf
    (
        fvc::interpolate
        (
            scalar(1.0)/alphaf,
            "interpolate(oneoveralphaf)"
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
        volVectorField HsbyA(constrainHbyA(rUsA*UsEqn.H(), Us, p_rgh));
        volVectorField HfbyA(constrainHbyA(rUfA*UfEqn.H(), Uf, p_rgh));

        surfaceScalarField phiHsbyA
        (
            "phiHsbyA",
            fvc::flux(HsbyA)
          + fvc::interpolate(rhos*rUsA)*fvc::ddtCorr(Us, phis)
        );

        surfaceScalarField phiHfbyA
        (
            "phiHfbyA",
            fvc::flux(HfbyA)
          + fvc::interpolate(rhof*rUfA)*fvc::ddtCorr(Uf, phif)
        );


        //adjustPhi(phiHsbyA, Us, p_rgh);

        // Momentum flux at the cell faces for phase a
        surfaceScalarField phiDrags
        (
          // grad(pa)
          - (oneOverAlphas*fvc::snGrad(ps))*rUsAf*mesh.magSf()
          // diagonal part of stress --> move to Rcs in UEqns.H ?
          + (oneOverAlphas*fvc::snGrad(lambdaUs*tr(gradUsT)))
           *rUsAf*mesh.magSf()
          // gravity
          + rUsAf
           *(
               fvc::interpolate(rhos-rhof)*(gravity & mesh.Sf())
            )
          // mean pressure gradient
          + rUsAf*(gradPOSC & mesh.Sf())
        );
        if (not SUSaniso)
        {
            //turbulent suspension term
            phiDrags += -(fvc::interpolate(SUS*K*turbulencef->nut()
                                *(alphaf/(alphas + alphasSmall)))
                                *fvc::snGrad(alphas))*rUsAf*mesh.magSf();
        }

        // Momentum flux at the cell faces for phase b
        surfaceScalarField phiDragf
        (
          // driving force
            rUfAf*(gradPOSC & mesh.Sf())
            //surface tension
          - (fvc::interpolate(sigma*Kg)*fvc::snGrad(alphag))*rUfAf*mesh.magSf()
          - (ghf*fvc::snGrad(rhof))*rUfAf*mesh.magSf()
          // Convective acceleration for 1D non-uniform simulation
          //- rUfAf*fvc::interpolate(neg(alphas-0.3))*((gradConvAcc/rhof)
          //  & mesh.Sf())
        );
        if (not SUSaniso)
        {
            // turbulent suspension term
            phiDragf += fvc::interpolate(SUS*K*
                        turbulencef->nut())*fvc::snGrad(alphas)
                        *rUfAf*mesh.magSf();
        }

        // Fix for gravity on outlet boundary.
        forAll(p_rgh.boundaryField(), patchi)
        {
             if (isA<zeroGradientFvPatchScalarField>
                 (p_rgh.boundaryField()[patchi]))
             {
                 phiDrags.boundaryFieldRef()[patchi] = 0.0;
                 phiDragf.boundaryFieldRef()[patchi] = 0.0;
             }
             if (isA<symmetryFvPatchScalarField>
                 (p_rgh.boundaryField()[patchi]))
             {
                 phiDrags.boundaryFieldRef()[patchi] = 0.0;
                 phiDragf.boundaryFieldRef()[patchi] = 0.0;
             }
             if (p_rgh.boundaryField().types()[patchi] ==
                 "symmetryPlane")
             {
                 phiDrags.boundaryFieldRef()[patchi] = 0.0;
                 phiDragf.boundaryFieldRef()[patchi] = 0.0;
             }
             if (p_rgh.boundaryField().types()[patchi] ==
                 "groovyBC")
             {
                 phiDrags.boundaryFieldRef()[patchi] = 0.0;
                 phiDragf.boundaryFieldRef()[patchi] = 0.0;
             }
             if (p_rgh.boundaryField().types()[patchi] ==
                 "outletInlet")
             {
                 phiDrags.boundaryFieldRef()[patchi] = 0.0;
                 phiDragf.boundaryFieldRef()[patchi] = 0.0;
             }
             if (p_rgh.boundaryField().types()[patchi] ==
                 "fixedFluxPressure")
             {
                 phiDrags.boundaryFieldRef()[patchi] = 0.0;
                 phiDragf.boundaryFieldRef()[patchi] = 0.0;
             }
        }

        phiHsbyA += phiDrags;
        phiHfbyA += phiDragf;

        // Face-drag coefficients
        surfaceScalarField rAUKds(fvc::interpolate(alphaf*K*rUsA));
        surfaceScalarField rAUKdf(fvc::interpolate(alphas*K*rUfA));
        // Compute predicted velocity fluxes at cell faces including
        //       including semi-implicit drag contribution
        surfaceScalarField phiHbyA
        (
            alphasf*(phiHsbyA + rAUKds*phif)
          + alphaff*(phiHfbyA + rAUKdf*phis)
        );
        // adjustPhi(phiHbyA, U, p_rgh);

        // Assemble Poisson operator for the pressure
        surfaceScalarField Dp
        (
            "(rho*(1|A(U)))",
            alphasf*rUsAf + alphaff*rUfAf
        );

        // Update the fixedFluxPressure BCs to ensure flux consistency
        setSnGrad<fixedFluxPressureFvPatchScalarField>
        (
            p_rgh.boundaryFieldRef(),
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
            fvScalarMatrix p_rghEqn
            (
                fvm::laplacian(Dp, p_rgh) == fvc::div(phiHbyA)
            );

            p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));
#if (defined(OFVERSION) && (OFVERSION >= 700))
            p_rghEqn.solve();
#else
            p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));
#endif
        if (pimple.finalNonOrthogonalIter())
        {
            // Compute pressure gradient at cell faces
                surfaceScalarField SfGradp(p_rghEqn.flux()/Dp);
                // Compute corrected velocity fluxes at cell faces
                //       due to pressure gradient
                surfaceScalarField phiss(phiHsbyA - rUsAf*SfGradp);
                surfaceScalarField phifs(phiHfbyA - rUfAf*SfGradp);
                // Compute mixture velocity flux at cell faces
                phi = phiHbyA - p_rghEqn.flux();
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

                p_rgh.relax();
                // Is it useful?
                SfGradp = p_rghEqn.flux()/Dp;

                // Velocity correction step at cell centers
                volVectorField Ufs
                (
                    HfbyA + rUfA*fvc::reconstruct(phiDragf/rUfAf-SfGradp)
                );
                volVectorField Uss
                (
                    HsbyA + rUsA*fvc::reconstruct(phiDrags/rUsAf-SfGradp)
                );
                // Compute corrected relative velocity due to drag force
                //       at cell centers
                volScalarField Ds(alphaf*K*rUsA);
                volScalarField Df(alphas*K*rUfA);
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
        }
    }

    #include "continuityErrs.H"

    p == p_rgh + rhof*gh;

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rhof*gh;
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
