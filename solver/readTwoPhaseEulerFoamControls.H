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
#include "readTimeControls.H"

int nGammaCorr(pimple.dict().get<int>("nGammaCorr"));
int nAlphasCorr(pimple.dict().get<int>("nAlphasCorr"));
int nAlphaSubCycles(pimple.dict().get<int>("nAlphaSubCycles"));
bool MULESCorr(pimple.dict().lookupOrDefault<Switch>("MULESCorr", false));
bool gammaApplyPrevCorr(pimple.dict().lookupOrDefault<Switch>("gammaApplyPrevCorr", false));
// compression coefficients
scalar cAlpha
(
	pimple.dict().lookupOrDefault<scalar>("cAlpha", 0)
);
scalar icAlpha
(
	pimple.dict().lookupOrDefault<scalar>("icAlpha", 0)
);
	 
// Shear compression coefficient
scalar scAlpha
(
	pimple.dict().lookupOrDefault<scalar>("scAlpha", 0)
);