{
    word scheme("div(phi,alphas)");
    word schemer("div(phir,alphas)");
    word gammaScheme("div(phi,gamma)");
    word gammarScheme("div(phir,gamma)");

    #include "updateSurfaceTension.H"

    surfaceScalarField phir("phir", phis - phif);
    // Standard face-flux compression coefficient
    surfaceScalarField phic(cAlpha*mag(phi/mesh.magSf()));

    //#include "updateSurfaceTension.H"
    for (int acorr=0; acorr<nAlphasCorr; acorr++)
    {
        // Create the limiter to be used for all phase-fractions
        scalarField allLambda(mesh.nFaces(), 1.0);

        // Split the limiter into a surfaceScalarField
        slicedSurfaceScalarField lambda
        (
            IOobject
            (
                "lambda",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimless,
            allLambda,
            false   // Use slices for the couples
        );

        surfaceScalarField phirc(phic*nHatf);

        surfaceScalarField gammaPhi
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

        // Create the bounded (upwind) flux for gamma
        surfaceScalarField gammaPhiBD
        (
            upwind<scalar>(mesh, phi).flux(gamma)
        );

        // Calculate the flux correction for gamma
        gammaPhi -= gammaPhiBD;

        MULES::limiter
        (
            allLambda,
            1.0/runTime.deltaT().value(),
            geometricOneField(),
            gamma,
            gammaPhiBD,
            gammaPhi,
            zeroField(),
            zeroField(),
            oneField(),
            zeroField()
        );

        gammaPhi = gammaPhiBD + lambda*gammaPhi;

        // Reset allLambda to 1.0
        allLambda = 1.0;

        // Solve for gamma
        solve(fvm::ddt(gamma) + fvc::div(gammaPhi));// - gamma*fvc::div(Uf));

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
}
