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
#include "createTwoPhaseTurbulence.H"

singlePhaseTransportModel phasesTransport(Us, phis);
volScalarField alphasEff
(
    IOobject
    (
        IOobject::groupName("alphasEff", "sol"),
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimless, Zero)
);

alphasEff = alphas+10*alphasSmall;
surfaceScalarField alphascf("alphascf", fvc::interpolate(alphasEff));
surfaceScalarField alphasPhis
(
    IOobject::groupName
    (
        "alphasPhis",
        "sol"
    ),
    alphascf*phis
);

autoPtr<sedIncompressibleTurbulenceModel<singlePhaseTransportModel>>
turbulences
(
    sedIncompressibleTurbulenceModel<singlePhaseTransportModel>::New
    (
        alphasEff,
        Us,
        alphasPhis,
        phis,
        phasesTransport
    )
);

singlePhaseTransportModel phasefTransport(Uf, phif);

surfaceScalarField alphafcf("alphafcf", fvc::interpolate(alphaf));
surfaceScalarField alphafPhif
(
    IOobject::groupName
    (
        "alphafPhi",
        "flu"
    ),
    alphafcf*phif
);

autoPtr<sedIncompressibleTurbulenceModel<singlePhaseTransportModel>>
turbulencef
(
    sedIncompressibleTurbulenceModel<singlePhaseTransportModel>::New
    (
        alphaf,
        Uf,
        alphafPhif,
        phif,
        phasefTransport
    )
);

IOdictionary turbulencePropertiess
(
    IOobject
    (
        "turbulenceProperties.sol",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

IOdictionary turbulencePropertiesf
(
    IOobject
    (
        "turbulenceProperties.flu",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

Info<< "Calculating field nuEffs\n" << endl;
volScalarField nuEffs
(
    IOobject
    (
        "nuEffs",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0, 2, -1, 0, 0), 0.0)
);
Info<< "Calculating field nuFrs\n" << endl;
volScalarField nuFrs
(
    IOobject
    (
        "nuFrs",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0, 2, -1, 0, 0), 0.0)
);

Info<< "Calculating field nuEfff\n" << endl;
volScalarField nuEfff
(
    IOobject
    (
        "nuEfff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0, 2, -1, 0, 0), 0.0)
);

volScalarField spherSigmaSGSf
(
    IOobject
    (
        "sphereSigmaSGSf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0, 2, -2, 0, 0), 0.0)
);

volScalarField spherSigmaSGSs
(
    IOobject
    (
        "spherSigmaSGSs",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0, 2, -2, 0, 0), 0.0)
);
