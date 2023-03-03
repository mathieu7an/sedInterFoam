/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Application
    sedInterFoam

Description
    Solver for a system of 3 phases with one phase dispersed,
    e.g.  solid particles in water and air.

Version
    1.0

Author
    Julien Chauchat, Cyrille Bonamy, Antoine Mathieu, Rémi Chassagne,
    Tim Nagel, Zhen Cheng, Tian-Jian Hsu and Eduard Puig Montella.

Date
    January 11, 2023

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "sedIncompressibleTurbulenceModel.H"
#include "unitConversion.H"

#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"

#include "symmetryFvPatchFields.H"
#include "fixedFluxPressureFvPatchScalarField.H"
#include "slipFvPatchFields.H"
#include "partialSlipFvPatchFields.H"
#include "alphaContactAngleTwoPhaseFvPatchScalarField.H"
#include "relaxationZone.H"

#include "dragModel.H"
#include "phaseModel.H"
#include "ppModel.H"

#include "kineticTheoryModel.H"
#include "granularRheologyModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "fixedValueFvsPatchFields.H"
//#include "IOMRFZoneList.H"
//#include "IOMRFZoneList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
 //   #include "postProcess.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"

    #include "readGravity.H"
    #include "readWaveProperties.H"
    #include "createFields.H"
    #include "createTurbulence.H"
    #include "createFvOptions.H"

    #include "initContinuityErrs.H"
    #include "createTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"
    #include "createFavreAveraging.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    // Test on SUSlocal
    //
    if (SUSlocal)
    {
        Info<< "\nLocal Schmidt number activated" << endl;
    }
    else
    {
       if (max(SUS).value() == 0)
       {
           Info<< "Turbulence suspension term is neglected" << endl;
       }
       else if (max(SUS).value() > 0)
       {
           Info<< "Turbulence suspension term is included" << endl;
       }
       else
       {
           Info<< "Turbulence suspension coefficient SUS can't be negative"
               << endl;
       }
    }
    // Test on granular stress model
    if (kineticTheory.on() && granularRheology.on())
    {
        Info<< "\nKinetic theory and granular rheology are set on." << endl;
        Info<< " This option is not supported!" << endl;
    }

    // stress formulation
    Switch faceMomentum
    (
        pimple.dict().lookupOrDefault<Switch>("faceMomentum", false)
    );
    Info<< "Choice for faceMomentum : "<<faceMomentum
        << endl;
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTwoPhaseEulerFoamControls.H"
        #include "CourantNos.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"

        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

//      Apply a ramp in time on the gravity acceleration
        #include "gravityRamp.H"

//      Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
//          Solve for mass conservation equations
            #include "gammaEqn.H"

            #include "alphaEqn.H"

            //#include "updateSurfaceTension.H"

//          Compute lift and drag coefficients
            #include "liftDragCoeffs.H"

//          Compute the granular stress: pff, nuFrs, nuEffs and lambdaUs
//             from Kinetic Theory of granular flows or mu(I) rheology
            #include "callGranularStress.H"

//          Assemble the momentum balance equations for both phases a and b
//          And assemble and solve the pressure poisson equation
//             and apply the velocity correction step for both phases a and b
            if (faceMomentum)
            {
                #include "pU/UEqns.H"
                #include "pU/pEqn.H"
            }
            else
            {
                #include "pU/UEqns.H"
                #include "pU/pEqn.H"
            }

//          Compute the phase accelerations for added mass force
            //#include "DDtU.H"

            if (pimple.turbCorr())
            {
//              Solve for turbulence models
                #include "updateTwoPhaseTurbulence.H"
                turbulencef->correct();
                if (turbulencePropertiesf.get<word>("simulationType")=="LES")
                {
                    spherSigmaSGSf = turbulencef->spherSigmaSGS();
                }
                turbulences->correct();
                if (turbulencePropertiess.get<word>("simulationType")=="LES")
                {
                    spherSigmaSGSs = turbulences->spherSigmaSGS();
                }

                if (debugInfo)
                {
                    Info << " max(nutf) = "
                         << max(turbulencef->nut()).value() << endl;
                }
            }
        }
        if (debugInfo)
        {
            Info<< "min(Us) = " << gMin(Us)
                << "max(Us) = " << gMax(Us) << endl;
            Info<< "min(Uf) = " << gMin(Uf)
                << "max(Uf) = " << gMax(Uf) << nl << endl;
        }
//      Write output
        #include "OutputGradPOSC.H"
        #include "writeOutput.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}
// ************************************************************************* //
