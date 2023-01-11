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
Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

IOdictionary twophaseRASProperties
(
    IOobject
    (
        "twophaseRASProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

IOdictionary ppProperties
(
    IOobject
    (
        "ppProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// Read ppProperties


dimensionedScalar alphasMax
(
    dimensionedScalar::getOrDefault
    (
        "alphasMax",
        ppProperties,
        dimless,
        0.635
    )
);

dimensionedScalar alphasMinFriction
(
    dimensionedScalar::getOrDefault
    (
        "alphasMinFriction",
        ppProperties,
        dimless,
        0.57
    )
);

dimensionedScalar Fr_
(
    dimensionedScalar::getOrDefault
    (
        "Fr",
        ppProperties,
        dimensionSet(1, -1, -2, 0, 0, 0, 0),
        5e-2
    )
);
dimensionedScalar eta0_
(
    dimensionedScalar::getOrDefault
    (
        "eta0",
        ppProperties,
        dimless,
        3
    )
);
dimensionedScalar eta1_
(
    dimensionedScalar::getOrDefault
    (
        "eta1",
        ppProperties,
        dimless,
        5
    )
);
Switch packingLimiter(ppProperties.get<Switch>("packingLimiter"));

// Read transportProperties

dimensionedScalar alphaDiffusion
(
    dimensionedScalar::getOrDefault
    (
        "alphaDiffusion",
        transportProperties,
        dimensionSet(0, 2, -1, 0, 0, 0, 0),
        0e0
    )
);

dimensionedScalar alphasSmall
(
    dimensionedScalar::getOrDefault
    (
        "alphasSmall",
        transportProperties,
        dimless,
        1e-6
    )
);

dimensionedScalar nuMax
(
    dimensionedScalar::getOrDefault
    (
        "nuMax",
        transportProperties,
        dimensionSet(0, 2, -1, 0, 0, 0, 0),
        1e0
    )
);

dimensionedScalar Dsmall
(
    dimensionedScalar::getOrDefault
    (
        "Dsmall",
        transportProperties,
        dimensionSet(0, 0, -1, 0, 0, 0, 0),
        1e-6
    )
);

dimensionedScalar nuMaxExp
(
    dimensionedScalar::getOrDefault
    (
        "nuMaxExp",
        transportProperties,
        dimensionSet(0, 2, -1, 0, 0, 0, 0),
        1e0
    )
);
// Initialize phase a and b variables

autoPtr<phaseModel> phases = phaseModel::New
(
    mesh,
    transportProperties,
    "sol"
);

autoPtr<phaseModel> phasew = phaseModel::New
(
    mesh,
    transportProperties,
    "wat"
);

autoPtr<phaseModel> phaseg = phaseModel::New
(
    mesh,
    transportProperties,
    "gas"
);

volVectorField& Us = phases->U();
volScalarField& alphas = phases->alpha();
surfaceScalarField& phis = phases->phi();
const dimensionedScalar& rhos = phases->rho();
const dimensionedScalar& nus = phases->nu();

volVectorField& Uw = phasew->U();
volScalarField& alphaw = phasew->alpha();
surfaceScalarField& phiw = phasew->phi();
const dimensionedScalar& rhow = phasew->rho();
const dimensionedScalar& nuw = phasew->nu();

volVectorField& Ug = phaseg->U();
volScalarField& alphag = phaseg->alpha();
surfaceScalarField& phig = phaseg->phi();
const dimensionedScalar& rhog = phaseg->rho();
const dimensionedScalar& nug = phaseg->nu();

volScalarField alphaf
(
    IOobject
    (
         IOobject::groupName("alpha", "flu"),
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0)
);

alphaf = 1.0 - alphas;

volScalarField gamma
(
    IOobject
    (
        "gamma",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

alphaw = gamma*alphaf;
alphag = 1.0 - alphaw - alphas;

volVectorField Uf
(
    IOobject
    (
        IOobject::groupName("U", "flu"),
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    Uw
);

Uf = (alphag*Ug + alphaw*Uw)/alphaf;

surfaceScalarField phif
(
    IOobject
    (
        IOobject::groupName("phi", "flu"),
        runTime.timeName(),
        mesh
    ),
    fvc::flux(Uf)
);

volScalarField rhof
(
    IOobject
    (
        "rhof",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::AUTO_WRITE
    ),
    (alphag*rhog + alphaw*rhow)/alphaf
);

// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alphaw*rhow + alphag*rhog + alphas*rhos
);
rho.oldTime();

volScalarField muf
(
    IOobject
    (
        "muf",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::AUTO_WRITE
    ),
    (alphaw*rhow*nuw + alphag*rhog*nug)
);

volScalarField nuf
(
    IOobject
    (
        "nuf",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::AUTO_WRITE
    ),
    muf/(rhof*alphaf)
);

// Reading ForceProperties
Info<< "\nReading forceProperties for gradP\n" <<endl;

IOdictionary forceProperties
(
    IOobject
    (
        "forceProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

dimensionedScalar tilt
(
    dimensionedScalar::getOrDefault
    (
        "tilt",
        forceProperties,
        dimless,
        0
    )
);

dimensionedVector gradPAMP1
(
    dimensionedVector::getOrDefault
    (
        "gradPAMP1",
        forceProperties,
        dimensionSet(1, -2, -2, 0, 0, 0, 0),
        vector(0, 0, 0)
    )
);

dimensionedVector gradPAMP2
(
    dimensionedVector::getOrDefault
    (
        "gradPAMP2",
        forceProperties,
        dimensionSet(1, -2, -2, 0, 0, 0, 0),
        vector(0, 0, 0)
    )
);

dimensionedVector gradPAMP3
(
    dimensionedVector::getOrDefault
    (
        "gradPAMP3",
        forceProperties,
        dimensionSet(1, -2, -2, 0, 0, 0, 0),
        vector(0, 0, 0)
    )
);

dimensionedVector gradPAMP4
(
    dimensionedVector::getOrDefault
    (
        "gradPAMP4",
        forceProperties,
        dimensionSet(1, -2, -2, 0, 0, 0, 0),
        vector(0, 0, 0)
    )
);

dimensionedVector gradPAMP5
(
    dimensionedVector::getOrDefault
    (
        "gradPAMP5",
        forceProperties,
        dimensionSet(1, -2, -2, 0, 0, 0, 0),
        vector(0, 0, 0)
    )
);

dimensionedVector gradPMEAN
(
    dimensionedVector::getOrDefault
    (
        "gradPMEAN",
        forceProperties,
        dimensionSet(1, -2, -2, 0, 0, 0, 0),
        vector(0, 0, 0)
    )
);

dimensionedVector gradPOSC
(
    dimensionedVector::getOrDefault
    (
        "gradPOSC",
        forceProperties,
        dimensionSet(1, -2, -2, 0, 0, 0, 0),
        vector(0, 0, 0)
    )
);

dimensionedVector gradConvAcc
(
    dimensionedVector::getOrDefault
    (
        "gradConvAcc",
        forceProperties,
        dimensionSet(1, -2, -2, 0, 0, 0, 0),
        vector(0, 0, 0)
    )
);
dimensionedScalar gradPOSCMin ("gradPOSCMin", gradPOSC.dimensions(), 1e-6);

dimensionedScalar oscpT
(
    dimensionedScalar::getOrDefault
    (
        "oscpT",
        forceProperties,
        dimless,
        0
    )
);

// to avoid division by 0 in gradP impose 1 if 0 is specified
dimensionedScalar oscpTMin ("oscpTMin", oscpT.dimensions(), 1e-6);
oscpT=max(oscpTMin, oscpT);

dimensionedScalar Cvm
(
    dimensionedScalar::getOrDefault
    (
        "Cvm",
        forceProperties,
        dimless,
        0
    )
);

dimensionedScalar Cl
(
    dimensionedScalar::getOrDefault
    (
        "Cl",
        forceProperties,
        dimless,
        0
    )
);

dimensionedScalar Ct
(
    dimensionedScalar::getOrDefault
    (
        "Ct",
        forceProperties,
        dimless,
        0
    )
);

dimensionedScalar initTheta1
(
    dimensionedScalar::getOrDefault
    (
        "initTheta1",
        forceProperties,
        dimless,
        0
    )
);

dimensionedScalar initTheta2
(
    dimensionedScalar::getOrDefault
    (
        "initTheta2",
        forceProperties,
        dimless,
        0
    )
);

dimensionedScalar initTheta3
(
    dimensionedScalar::getOrDefault
    (
        "initTheta3",
        forceProperties,
        dimless,
        0
    )
);

dimensionedScalar initTheta4
(
    dimensionedScalar::getOrDefault
    (
        "initTheta4",
        forceProperties,
        dimless,
        0
    )
);

dimensionedScalar initTheta5
(
    dimensionedScalar::getOrDefault
    (
        "initTheta5",
        forceProperties,
        dimless,
        0
    )
);

dimensionedScalar InitFreeze
(
    dimensionedScalar::getOrDefault
    (
        "InitFreeze",
        forceProperties,
        dimless,
        0
    )
);

dimensionedScalar fixbeddepth
(
    dimensionedScalar::getOrDefault
    (
        "fixbeddepth",
        forceProperties,
        dimensionSet(0, 1, 0, 0, 0, 0, 0),
        0
    )
);

dimensionedScalar ClipUs
(
    dimensionedScalar::getOrDefault
    (
        "ClipUs",
        forceProperties,
        dimless,
        0
    )
);

// Debuging info
Switch debugInfo
(
    forceProperties.getOrDefault<Switch>("debugInfo", false)
);
// Sponge Layer
Switch spongeLayer
(
    forceProperties.getOrDefault<Switch>("spongeLayer", false)
);
dimensionedScalar xSmin
(
    dimensionedScalar::getOrDefault
    (
        "xSmin",
        forceProperties,
        dimensionSet(0, 1, 0, 0, 0, 0, 0),
        0
    )
);
dimensionedScalar xSmax
(
    dimensionedScalar::getOrDefault
    (
        "xSmax",
        forceProperties,
        dimensionSet(0, 1, 0, 0, 0, 0, 0),
        0
    )
);

Switch writeTau
(
    forceProperties.getOrDefault<Switch>("writeTau", false)
);


Info<< "Initializing with specified pressure gradient:"
<< gradPOSC.value() << endl << endl;

// open files to write gradP info
OFstream gradPOSCFile
(
    "gradPOSC.txt"
);

// Initialize variables
volScalarField Unity
(
    IOobject
    (
        "Unity",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("one", dimensionSet(0, 0, 0, 0, 0), 1.0)
);

Info<< "Calculating field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    alphas*Us + alphaf*Uf + alphag*Ug
);

Info<< "Calculating face flux field phi\n" << endl;
surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::flux(U)
);

#include "UsUfLimiter.H"
Us.correctBoundaryConditions();

// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);

Info<< "Calculating field DDtUs and DDtUf\n" << endl;

volVectorField DDtUs
(
    fvc::ddt(Us)
  + fvc::div(phis, Us)
  - fvc::div(phis)*Us
);

volVectorField DDtUf
(
    fvc::ddt(Uf)
  + fvc::div(phif, Uf)
  - fvc::div(phif)*Uf
);

// MULES flux from previous time-step
surfaceScalarField gammaPhi
(
    IOobject
    (
        "gammaPhi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    phi*fvc::interpolate(gamma)
);

surfaceScalarField alphasPhi
(
    IOobject
    (
        "alphasPhi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    phi*fvc::interpolate(alphas)
);


dimensionedScalar sigmagf
(
    dimensionedScalar::getOrDefault
    (
        "sigmagf",
        transportProperties,
        dimensionSet(1, 0, -2, 0, 0, 0, 0),
        0.05
    )
);
dimensionedScalar sigmags
(
    dimensionedScalar::getOrDefault
    (
        "sigmags",
        transportProperties,
        dimensionSet(1, 0, -2, 0, 0, 0, 0),
        0.04
    )
);

// MULES Correction
tmp<surfaceScalarField> tgammaPhiCorr0;

//The following two flags are used to solve the stability problem which occurs when air is entrained
//on water-sediment interface. It may cause mass conservation error which can be large in some cases 

//Flag to remove air from water-sand interface
Switch removeAir(transportProperties.lookupOrDefault("removeAir", false));

//Flag to add water instead of the removed air
Switch addWater(transportProperties.lookupOrDefault("addWater", false)); 

volScalarField sigma
(
    IOobject
    (
        "sigma",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(1, 0, -2, 0, 0), 0.0)
);

volScalarField Kg
(
    IOobject
    (
        "Kg",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0, -1, 0, 0, 0), 0.0)
);

// Cell gradient of alphaw
volVectorField gradAlphaw(fvc::grad(alphaw, "nHat"));;

// Interpolated face-gradient of alphaw
surfaceVectorField gradAlphawf(fvc::interpolate(gradAlphaw));

dimensionedScalar deltaN(1e-8/cbrt(average(mesh.V())));

surfaceScalarField nHatf
(
    IOobject
    (
        "nHatf",
        runTime.timeName(),
        mesh
    ),
    (gradAlphawf/(mag(gradAlphawf) + deltaN)) & mesh.Sf()
);

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
Info<< "Reading field ps\n" << endl;
volScalarField ps
(
    IOobject
    (
        "ps",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(1, -1, -2, 0, 0), 0.0)
);

volTensorField gradUsT
(
    fvc::grad(Us)().T()
);

volTensorField gradUfT
(
    fvc::grad(Uf)().T()
);

volTensorField Tauw
(
    IOobject
    (
        "Tauw",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rhow*(scalar(1.0)-alphas)*nuw*dev(fvc::grad(Uw))
);

volTensorField Taus
(
    IOobject
    (
        "Taus",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    rhos*alphas*nus*dev(fvc::grad(Us))
);

volTensorField Taug
(
    IOobject
    (
        "Taug",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    rhog*alphag*nug*dev(fvc::grad(Ug))
);

volTensorField Tauf
(
    IOobject
    (
        "Tauf",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    rhof*alphaf*nuf*dev(fvc::grad(Uf))
);

//creating alphasPlastic field
Info<< "Computing field alphasPlastic\n" << endl;
volScalarField alphasPlastic
(
     IOobject
     (
         "alphasPlastic",
         runTime.timeName(),
         mesh,
         IOobject::READ_IF_PRESENT,
         IOobject::NO_WRITE
     ),
     mesh,
     dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0)
);

volScalarField delta
(
    IOobject
    (
        "delta",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0)
);

Info<< "Calculating field g.h\n" << endl;
volScalarField gh("gh", gravity & mesh.C());
surfaceScalarField ghf("ghf", gravity & mesh.Cf());

// Read interfacial properties

IOdictionary interfacialProperties
(
    IOobject
    (
        "interfacialProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

autoPtr<dragModel> drags = dragModel::New
(
    interfacialProperties,
    phases,
    rhof,
    nuf
);

autoPtr<dragModel> dragf = dragModel::New
(
    interfacialProperties,
    phases,
    rhof,
    nuf
);

word dragPhase("blended");
if (interfacialProperties.found("dragPhase"))
{
    dragPhase = interfacialProperties.get<word>("dragPhase");

    bool validDrag =
        dragPhase == "sol" || dragPhase == "flu" || dragPhase == "blended";

    if (validDrag == false)
    {
        FatalErrorIn(args.executable())
            << "invalid dragPhase " << dragPhase
            << exit(FatalError);
    }
}


Info << "dragPhase is " << dragPhase << endl;

// Compute fall velocity
#include "calculateUfall.H"
dimensionedScalar Ufall
(
   "Ufall",
   dimensionSet(0, 1, -1, 0, 0, 0, 0),
   max(UfallField).value()
);

// Initalize particle pressure model
autoPtr<ppModel> pp = ppModel::New
(
    ppProperties,
    phases,
    phasew
);

// Initalize Kinetic theory
kineticTheoryModel kineticTheory
(
    phases,
    Uf,
    drags
);

// Initialize granular rheology
granularRheologyModel granularRheology
(
    phases,
    rhof,
    nuf,
    ps,
    Dsmall
);

volScalarField lambdaUs
(
    IOobject
    (
        "lambdaUs",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(1, -1, -1, 0, 0), 0.0)
);

volScalarField K
(
    IOobject
    (
        "K",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    drags->K(mag(Uf-Us))
);

volVectorField gradPs
(
    IOobject
    (
        "gradPs",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimensionSet(0, 1, -2, 0, 0), vector::zero)
);

volVectorField gradPf
(
    IOobject
    (
        "gradPf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimensionSet(0, 1, -2, 0, 0), vector::zero)
);

volVectorField suss
(
    IOobject
    (
        "suss",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimensionSet(0, 1, -2, 0, 0), vector::zero)
);
volVectorField susf
(
    IOobject
    (
        "susf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimensionSet(0, 1, -2, 0, 0), vector::zero)
);
volScalarField pff
(
    IOobject
    (
        "pff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(1, -1, -2, 0, 0), 0.0)
);

volScalarField SUStilde
(
    IOobject
    (
        "SUStilde",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("one", dimensionSet(0, 0, 0, 0, 0), 1.0)
);

volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    p_rgh + rhof*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    mesh.solutionDict().subDict("PIMPLE"),
    pRefCell,
    pRefValue
);
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
Info<<"end of createFields"<<endl;
