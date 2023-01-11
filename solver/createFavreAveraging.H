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

Switch favreAveraging
(
    runTime.controlDict().lookupOrDefault("favreAveraging", false)
);

Switch TKEBudgetLES
(
    runTime.controlDict().lookupOrDefault("TKEBudgetLES", false)
);

scalar avCount = 0;

volVectorField UwPrimeF
(
    IOobject
    (
        "UwPrimeF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "UwPrimeF",
        dimVelocity,
        vector (0, 0, 0)
    )
);

volVectorField UsPrimeF
(
    IOobject
    (
        "UsPrimeF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "UsPrimeF",
        dimVelocity,
        vector (0, 0, 0)
    )
);

volScalarField pMeanF
(
    IOobject
    (
        "pMeanF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "pMeanF",
        dimPressure,
        0
    )
);

volScalarField pPrimeF
(
    IOobject
    (
        "pPrimeF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "pPrimeF",
        dimPressure,
        0
    )
);

volScalarField RS2
    (
        IOobject
        (
            "RS2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "RS2",
            dimensionSet(0, 2, -2, 0, 0, 0, 0),
            0.0
        )
    );

volScalarField RS4
    (
        IOobject
        (
            "RS4",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "RS4",
            dimensionSet(0, 2, -2, 0, 0, 0, 0),
            0.0
        )
    );

volScalarField RS2Mean
    (
        IOobject
        (
            "RS2Mean",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "RS2Mean",
            dimensionSet(0, 2, -2, 0, 0, 0, 0),
            0.0
        )
    );

volScalarField RS4Mean
    (
        IOobject
        (
            "RS4Mean",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "RS4Mean",
            dimensionSet(0, 2, -2, 0, 0, 0, 0),
            0.0
        )
    );

volTensorField alphasSqrUsMean
(
    IOobject
    (
        "alphasSqrUsMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedTensor
    (
        "alphasSqrUsMean",
        dimVelocity*dimVelocity,
        tensor::zero
    )
);

volTensorField alphawSqrUwMean
(
    IOobject
    (
        "alphawSqrUwMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedTensor
    (
        "alphawSqrUwMean",
        dimVelocity*dimVelocity,
        tensor::zero
    )
);

volVectorField alphawUwMean
    (
        IOobject
        (
            "alphawUwMean",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector
        (
            "alphawUwMean",
            dimVelocity,
            vector (0, 0, 0)
        )
    );

volVectorField alphasUsMean
    (
        IOobject
        (
            "alphasUsMean",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector
        (
            "alphasUsMean",
            dimVelocity,
            vector (0, 0, 0)
        )
    );

volScalarField alpha_sMean
    (
        IOobject
        (
            "alpha_sMean",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "alphasUsMean",
            dimensionSet(0, 0, 0, 0, 0, 0, 0),
            0.0
        )
    );

volVectorField UwMeanF
(
    IOobject
    (
        "UwMeanF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "UwMeanF",
        dimVelocity,
        vector (0, 0, 0)
    )
);

volVectorField UsMeanF
(
    IOobject
    (
        "UsMeanF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "UsMeanF",
        dimVelocity,
        vector (0, 0, 0)
    )
);
volTensorField UwPrime2MeanF
(
    IOobject
    (
        "UwPrime2MeanF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedTensor
    (
        "UwPrime2MeanF",
        dimVelocity*dimVelocity,
        tensor::zero
    )
);

volTensorField UsPrime2MeanF
(
    IOobject
    (
        "UsPrime2MeanF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedTensor
    (
        "UsPrime2MeanF",
        dimVelocity*dimVelocity,
        tensor::zero
    )
);

volScalarField TKEMean
(
    IOobject
    (
        "TKEMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "TKEMean",
        dimVelocity*dimVelocity,
        0
    )
);
volScalarField TKEMeanProd
(
    IOobject
    (
        "TKEMeanProd",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "TKEMeanProd",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0
    )
);

volSymmTensorField SGSstrainTensor
(
    IOobject
    (
        "SGSstrainTensor",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedSymmTensor
    (
        "SGSstrainTensor",
        dimVelocity/dimLength,
        symmTensor::zero
    )
);
volSymmTensorField BB
(
    IOobject
    (
        "BB",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedSymmTensor
    (
        "BB",
        dimVelocity/dimLength,
        symmTensor::zero
    )
);
volScalarField turbDiffusionMean
(
    IOobject
    (
        "turbDiffusionMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "turbDiffusionMean",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0.0
    )
);
volVectorField turbDiffMean
(
    IOobject
    (
        "turbDiffMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "turbDiffMean",
        dimVelocity*dimVelocity*dimVelocity,
        vector::zero
    )
);
volScalarField pressDiffusionMean
(
    IOobject
    (
        "pressDiffusionMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "pressDiffusionMean",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0.0
    )
);

volSymmTensorField strainTensor
(
    IOobject
    (
        "strainTensor",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedSymmTensor
    (
        "strainTensor",
        dimVelocity/dimLength,
        symmTensor::zero
    )
);

volScalarField viscDissMean
(
    IOobject
    (
        "viscDissMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "viscDissMean",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0.0
    )
);
volScalarField viscDiffMean
(
    IOobject
    (
        "viscDiffMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "viscDiffMean",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0.0
    )
);
volSymmTensorField SGSStrainTensor
(
    IOobject
    (
        "SGSstrainTensor",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedSymmTensor
    (
        "SGSStrainTensor",
        dimVelocity/dimLength,
        symmTensor::zero
    )
);

volScalarField SGSDissMean
(
    IOobject
    (
        "SGSDissMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "SGSDissMean",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0.0
    )
);
volVectorField SGSDiffMean
(
    IOobject
    (
        "SGSDiffMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "SGSDiffMean",
        dimVelocity*dimVelocity*dimVelocity,
        vector::zero
    )
);
volScalarField SGSDiffusionMean
(
    IOobject
    (
        "SGSDiffusionMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "SGSDiffusionMean",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0.0
    )
);
volScalarField dragContwMean
(
    IOobject
    (
        "dragContwMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "dragContwMean",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0.0
    )
);

volVectorField alphawFluct
(
    IOobject
    (
        "alphawFluct",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "alphawFluct",
        dimVelocity*dimVelocity*dimVelocity,
        vector::zero
    )
);

volVectorField alphawFluctMean
(
    IOobject
    (
        "alphawFluctMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "alphawFluctMean",
        dimVelocity*dimVelocity*dimVelocity,
        vector::zero
    )
);

// Checks if field alphawFluctMean is present
IOobject alphawFluctMeanHeader
    (
        "alphawFluctMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );

if (alphawFluctMeanHeader.typeHeaderOk<volVectorField>(true))
{
    volVectorField alphawFluctMean2
        (
            IOobject
            (
                "alphawFluctMean",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );
    alphawFluctMean = alphawFluctMean2;
}

// Checks if field alpha_aMean is present
IOobject alpha_sMeanHeader
    (
        "alpha_sMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );

if (alpha_sMeanHeader.typeHeaderOk<volScalarField>(true))
{
    avCount = runTime.timeIndex();
    Info << avCount << endl;
    volScalarField alpha_sMean2
        (
            IOobject
            (
                "alpha_sMean",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );
    alpha_sMean = alpha_sMean2;
    // Checks if field RS2Mean is present
    IOobject RS2MeanHeader
    (
        "RS2Mean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (RS2MeanHeader.typeHeaderOk<volScalarField>(true))
    {
        volScalarField RS2Mean2
        (
            IOobject
            (
                "RS2Mean",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );
        RS2Mean = RS2Mean2;
    }
    // Checks if field RS4Mean is present
    IOobject RS4MeanHeader
    (
        "RS4Mean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (RS4MeanHeader.typeHeaderOk<volScalarField>(true))
    {
        volScalarField RS4Mean2
        (
            IOobject
            (
                "RS4Mean",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );
        RS4Mean = RS4Mean2;
    }
    // Checks if field UwMeanF is present
    IOobject UwMeanFHeader
    (
        "UwMeanF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (UwMeanFHeader.typeHeaderOk<volVectorField>(true))
    {
        volVectorField UwMeanF2
        (
            IOobject
            (
                "UwMeanF",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );
        UwMeanF = UwMeanF2;
        alphawUwMean = UwMeanF*(1-alpha_sMean);
    }

    // Checks if field UsMeanF is present
    IOobject UsMeanFHeader
    (
        "UsMeanF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (UsMeanFHeader.typeHeaderOk<volVectorField>(true))
    {
        volVectorField UsMeanF2
        (
            IOobject
            (
                "UsMeanF",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );
        UsMeanF = UsMeanF2;
        alphasUsMean = UsMeanF*(alpha_sMean+1e-20);
    }

    // Checks if field UwPrime2MeanF is present
    IOobject UwPrime2MeanFHeader
    (
        "UwPrime2MeanF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (UwPrime2MeanFHeader.typeHeaderOk<volTensorField>(true))
    {
        volTensorField UwPrime2MeanF2
        (
            IOobject
            (
                "UwPrime2MeanF",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );
        UwPrime2MeanF = UwPrime2MeanF2;
        alphawSqrUwMean = (UwPrime2MeanF+UwMeanF*UwMeanF)*(1-alpha_sMean);
    }

    // Checks if field UsPrime2MeanF is present
    IOobject UsPrime2MeanFHeader
    (
        "UsPrime2MeanF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (UsPrime2MeanFHeader.typeHeaderOk<volTensorField>(true))
    {
        volTensorField UsPrime2MeanF2
        (
            IOobject
            (
                "UsPrime2MeanF",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );
        UsPrime2MeanF = UsPrime2MeanF2;
        alphasSqrUsMean = (UsPrime2MeanF+UsMeanF*UsMeanF)
                         *(alpha_sMean+1e-20);
    }
}
