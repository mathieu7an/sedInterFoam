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
// if alpha < alphasSmall then set Us = Uw
if (ClipUs.value()==1)
{
    surfaceScalarField alphasf(fvc::interpolate(alphas));
    forAll(alphas, cell)
    {
        Usclip[cell] = Us[cell];
        if (alphas[cell]<alphasSmall.value())
        {
            Usclip[cell] = Uf[cell];
        }
    }
    forAll(alphasf, face)
    {
        phisclip[face] = phis[face];
        if (alphasf[face]<alphasSmall.value())
        {
            phisclip[face] = phif[face];
        }
    }
    //phi = alphasf*phis + (1-alphasf)*phif;
}

if (InitFreeze.value()==1)
{
// in this way the bed is immobile in streamwise
    volScalarField temp(Us & gravity);
// extract vertical component of velocity component,we get ( 0 ,0,w)
    volVectorField Ustemp(temp*gravity/magSqr(gravity));

    scalar indv1((vector(0, 1, 0) & gravity.value()));
    scalar indv2((vector(0, 0, 1) & gravity.value()));

    int ind = 0;
    if (indv1 != 0) {ind = 1;}
    if (indv2 != 0) {ind = 2;}

    forAll(alphas, cell)
    {
// initially, we need to fix the bed level, and let the flow develop
// only cell centers are considered
// 0-->x, 1-->y, 2-->z
        if (mesh.cellCentres()[cell][ind] < fixbeddepth.value())
        {
// exclude the horizontal conponents
            Us[cell] = Ustemp[cell];
        }
    }
    phis = (fvc::interpolate(Us) & mesh.Sf());
    phi = fvc::interpolate(alphas)*phis + fvc::interpolate(alphaf)*phif;
}

if (InitFreeze.value()==2)
{
// in this way the bed is immobile in streamwise
    volScalarField temp(Us & gravity);
// extract vertical component of velocity component,we get ( 0 ,0,w)
    volVectorField Ustemp(temp*gravity/magSqr(gravity));
// extract spanwise component of velocity, we get ( 0,v,0)
    volVectorField Ustemp2((Us & vector(0, 1, 0))*vector(0, 1, 0));

    scalar indv1((vector(0, 1, 0) & gravity.value()));
    scalar indv2((vector(0, 0, 1) & gravity.value()));

    int ind = 0;
    if (indv1 != 0) {ind = 1;}
    if (indv2 != 0) {ind = 2;}

    forAll(alphas, cell)
    {
// initially, we need to fix the bed level, and let the flow develop
// only cell centers are considered
// 0-->x, 1-->y, 2-->z
       if (mesh.cellCentres()[cell][ind] < fixbeddepth.value())
       {
// only exclude the streamwise conponents
           Us[cell] = Ustemp[cell]+Ustemp2[cell];
       }
    }
    phis = (fvc::interpolate(Us) & mesh.Sf());
    phi = fvc::interpolate(alphas)*phis + fvc::interpolate(alphaf)*phif;
}
