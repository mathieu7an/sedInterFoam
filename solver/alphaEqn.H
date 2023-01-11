{
    word scheme("div(phi,alphas)");
    word schemer("div(phir,alphas)");

    surfaceScalarField phir("phir", phis - phif);

    //#include "updateSurfaceTension.H"
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
		if(nGammaCorr>0)
		{
			alphaw = alphaf*gamma;
		    alphag = min(max(alphaf*(1.0-gamma), scalar(0)), scalar(1.0));
		}
		else
		{
			alphaw = alphaf;
			alphag = 0.0;
		}       

        if (debugInfo)
        {
            Info<< "Dispersed phase volume fraction = "
            << alphas.weightedAverage(mesh.V()).value()
            << "  Min(alphas) = " << gMin(alphas)
            << "  Max(alphas) = " << gMax(alphas)
            << endl;
        }
    }
}
