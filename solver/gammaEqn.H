{
    word gammaScheme("div(phi,gamma)");
    word gammarScheme("div(phir,gamma)");

    surfaceScalarField phir("phir", phis - phif);
    // Standard face-flux compression coefficient
    surfaceScalarField phic(cAlpha*mag(phi/mesh.magSf()));

    // Add the optional isotropic compression contribution
    if (icAlpha > 0)
    {
        phic *= (1.0 - icAlpha);
        phic += (cAlpha*icAlpha)*fvc::interpolate(mag(U));
    }
    surfaceScalarField::Boundary& phicBf =
        phic.boundaryFieldRef();

    // Do not compress interface at non-coupled boundary faces
    // (inlets, outlets etc.)
    forAll(phic.boundaryField(), patchi)
    {
        fvsPatchScalarField& phicp = phicBf[patchi];

        if (!phicp.coupled())
        {
            phicp == 0;
        }
    }
    //### Alphaf calaculations ####
	if(nGammaCorr>0) //set nGammaCorr= 0 for problems without free surface 
	{
        if (MULESCorr)
        {
            fvScalarMatrix gammaEqn
            (
                
              fv::EulerDdtScheme<scalar>(mesh).fvmDdt(gamma)
            + fv::gaussConvectionScheme<scalar>
                (
                    mesh,   
                    phi,
                    upwind<scalar>(mesh, phi)
                ).fvmDiv(phi, gamma)
            /*+ fv::gaussConvectionScheme<scalar>
                (
                    mesh,   
                    fvc::flux(-phir, alphas, gammarScheme),
                    upwind<scalar>(mesh, fvc::flux(-phir, alphas, gammarScheme))
                ).fvmDiv(fvc::flux(-phir, alphas, gammarScheme), gamma)*/
            //-gamma*fvc::div(Uf)
            //- fvm::Sp(fvc::div(Uf)/(1-alphas), gamma)
            );

            gammaEqn.solve();
            //gamma.min(1.0);
            //gamma.max(0);
            //alphaw = gamma*alphaf;

            Info<< "Carrier phase volume fraction = "
                << alphaw.weightedAverage(mesh.Vsc()).value()
                << "  Min(" << alphaw.name() << ") = " << min(alphaw).value()
                << "  Max(" << alphaw.name() << ") = " << max(alphaw).value()
                << endl;

            tmp<surfaceScalarField> tgammaPhiUD(gammaEqn.flux());
            gammaPhi = tgammaPhiUD();
            if (gammaApplyPrevCorr && tgammaPhiCorr0.valid())
            {
                Info<< "Applying the previous iteration compression flux" << endl;
                MULES::correct(geometricOneField(), gamma, gammaPhi, tgammaPhiCorr0.ref(), oneField(), zeroField());

                gammaPhi += tgammaPhiCorr0();
            }

            // Cache the upwind-flux
            tgammaPhiCorr0 = tgammaPhiUD;

            //alphaf = min(max(1.0-alphas, scalar(0)), scalar(1.0));
			//alphaw = alphaf*gamma;
			//alphag = min(max(alphaf*(1.0-gamma), scalar(0)), scalar(1.0));

            //Replace air in water-sand interface by water
            //if(addWater)
            //{
            //    alphaw += alphag*(pos(alphas-alphag)*pos(alphaw-0.1)*pos(0.99-alphaw));
            //}

            //Remove air from water-sand interface
            //if(removeAir)
            //{
            //    alphag *= (1-pos(alphas-alphag)*pos(alphaw-0.1)*pos(0.99-alphaw));
            //}
        }
        for (int aCorr=0; aCorr<nGammaCorr; aCorr++)
        {
            #include "updateSurfaceTension.H"

            surfaceScalarField phirc(phic*nHatf);

            tmp<surfaceScalarField> tgammaPhiUn
            (
                fvc::flux
                (
                    phi,
                    gamma,
                    gammaScheme
                )
            + fvc::flux
                (
                -fvc::flux(-phirc, (1.0-gamma), gammarScheme),
                    gamma,
                    gammarScheme
                )
            /*+ fvc::flux
                (
                fvc::flux(-phir, alphas, gammarScheme),
                    gamma,
                    gammarScheme
                )*/
            );

            if (MULESCorr)
            {
                tmp<surfaceScalarField> tgammaPhiCorr(tgammaPhiUn() - gammaPhi);
                volScalarField gamma0("gamma0", gamma);

                MULES::correct(geometricOneField(), gamma, tgammaPhiUn(), tgammaPhiCorr.ref(), oneField(), zeroField());

                // Under-relax the correction for all but the 1st corrector
                if (aCorr == 0)
                {
                    gammaPhi += tgammaPhiCorr();
                }
                else
                {
                    gamma = 0.5*gamma + 0.5*gamma0;
                    gammaPhi += 0.5*tgammaPhiCorr();
                }
            }
            else
            {
                gammaPhi = tgammaPhiUn;

                MULES::explicitSolve(geometricOneField(), gamma, phi, gammaPhi, oneField(), zeroField());
            }

            //solve(fvm::ddt(gamma) + fvc::div(gammaPhi) - fvm::Sp(fvc::div(Uf)/(1-alphas), gamma));

            //gamma.min(1.0);
            //gamma.max(0);
            
            alphaf = 1.0-alphas;
			alphaw = alphaf*gamma;
			alphag = alphaf*(1.0-gamma);
        }
        if (gammaApplyPrevCorr && MULESCorr)
        {
            tgammaPhiCorr0 = gammaPhi - tgammaPhiCorr0;
        }

        Info<< "Carrier phase volume fraction = "
            << alphaw.weightedAverage(mesh.Vsc()).value()
            << "  Min(" << alphaw.name() << ") = " << min(alphaw).value()
            << "  Max(" << alphaw.name() << ") = " << max(alphaw).value()
            << endl;
    }
    //For problems without free surface: 
	else
	{
		Info<< "Free surface calculation is turned OFF" << endl << endl;
        gamma = 1.0;
		alphaf = min(max(1.0-alphas, scalar(0)), scalar(1.0));
		alphaw = alphaf;
		alphag = 0.0;
	}
    // Construct the comp   lete mass flux
    rhoPhi =
        gammaPhi*(rhow-rhog)
        + phif*rhog 
        + fvc::flux(-phir, gamma*alphas, gammarScheme)*(rhow-rhog);

    muf = (alphaw*rhow*nuw + alphag*rhog*nug)/alphaf;
    nuf = (alphaw*nuw + alphag*nug)/alphaf;
    rhof = (alphaw*rhow + alphag*rhog)/alphaf;
    rho = (alphaw*rhow + alphag*rhog + alphas*rhos);

}
