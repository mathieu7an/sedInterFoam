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
    object      nut;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 2 -1 0 0 0 0];

internalField   uniform 0.; 

boundaryField
{
    leftWall           
    {
        type            zeroGradient;
    }

    rightWall
    {
	type		zeroGradient;
    }

    lowerWallSmooth           
    {
        type            nutUSpaldingWallFunction;
        value           uniform 0;
    }
    lowerWallRough           
    {
        type            nutUSpaldingWallFunction;
        value           uniform 0;
    }
    atmosphere      
    {
        type            zeroGradient;
    }

    cbackWall
    {
        type            zeroGradient;
    }
    cfrontWall
    {
        type            zeroGradient;
    }
}


// ************************************************************************* //
