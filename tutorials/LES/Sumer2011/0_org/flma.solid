/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2.4.0                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      flma;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 4 -4 0 0 0 0];

internalField   uniform 0.1;

boundaryField
{
    rightWall
    {
        type            fixedValue;
        value           uniform 0;
    }
    leftWall
    {
        type            fixedValue;
        value           uniform 0;
    }
    atmosphere
    {
        type            zeroGradient;
    }
    lowerWallSmooth
    {
        type            fixedValue;
        value           uniform 0;
    }
    lowerWallRough
    {
        type            fixedValue;
        value           uniform 0;
    }
    cfrontWall
    {
        type            fixedValue;
        value           uniform 0;
    }
    cbackWall
    {
        type            fixedValue;
        value           uniform 0;
    }
}


// ************************************************************************* //
