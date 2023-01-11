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
if (packingLimiter)
{
    //since the particle pressure near the alphasMax is extremely large
    //we need to limit the particle concentration to alphasMax-smallnumber
    // Calculating exceeding volume fractions
    volScalarField alphasEx(max(alphas - alphasMax+scalar(0.005), scalar(0)));

    // Finding neighbouring cells of the whole domain
    labelListList neighbour = mesh.cellCells();
    scalarField cellVolumes(mesh.cellVolumes());

    forAll(alphasEx, celli)
    {
        // Finding the labels of the neighbouring cells
        labelList neighbourCell = neighbour[celli];

        // Initializing neighbouring cells contribution
        scalar neighboursEx = 0.0;

        forAll(neighbourCell, cellj)
        {
            labelList neighboursNeighbour = neighbour[neighbourCell[cellj]];
            scalar neighboursNeighbourCellVolumes = 0.0;

            forAll(neighboursNeighbour, cellk)
            {
                neighboursNeighbourCellVolumes +=
                    cellVolumes[neighboursNeighbour[cellk]];
            }
            neighboursEx +=
                alphasEx[neighbourCell[cellj]]*cellVolumes[celli]
               /neighboursNeighbourCellVolumes;
        }
        alphas[celli] += neighboursEx - alphasEx[celli];
    }
    alphas.correctBoundaryConditions();
}
