{
    word scheme("div(phi,alphas)");
    word schemer("div(phir,alphas)");

    //#include "updateSurfaceTension.H"

    surfaceScalarField phir("phir", phis - phif);

    for (int acorr=0; acorr<nAlphasCorr; acorr++)
    {

        alphas.correctBoundaryConditions();
        fvScalarMatrix alphasEqn
        (
             fvm::ddt(alphas)
           + fvm::div(phi, alphas, scheme)
           + fvm::div(-fvc::flux(-phir, (1.0-alphas), schemer), alphas, schemer)
        );

        alphasEqn.relax();
        alphasEqn.solve();

        if (packingLimiter)
        {
            int i=0;
            while (i<100 && max(alphas).value()>=alphasMax.value())
            {
                #include "packingLimiter.H"
                i++;
            }
        }
        alphas.min(alphasMax);
        alphas.max(0);

		alphaf = min(max(1.0-alphas, scalar(0)), scalar(1.0));
        alphaw = alphaf*gamma;
        voidFrac = 1.0 - gamma;
		alphag = min(max(1.0 - alphaw - alphas, scalar(0)), scalar(1.0));

        //alphaw += alphag*(pos(alphas-alphag)*pos(alphaw-0.1)*pos(0.99-alphaw));
        //alphag *= (1-pos(alphas-alphag)*pos(alphaw-0.1)*pos(0.99-alphaw));

        if (debugInfo)
        {
            Info<< "Dispersed phase volume fraction = "
            << alphas.weightedAverage(mesh.V()).value()
            << "  Min(alphas) = " << gMin(alphas)
            << "  Max(alphas) = " << gMax(alphas)
            << endl;
        }
        // Construct the comp   lete mass flux
        /*rhoPhi =
            gammaPhi*(rhow-rhog)
            + phif*rhog
            + fvc::flux(-phir, gamma*alphas, gammarScheme)*(rhow-rhog);*/

        muf = (alphaw*rhow*nuw + alphag*rhog*nug)/alphaf;
        nuf = (alphaw*nuw + alphag*nug)/alphaf;
        rhof = (alphaw*rhow + alphag*rhog)/alphaf;
        rho = (alphaw*rhow + alphag*rhog + alphas*rhos);
    }
}
