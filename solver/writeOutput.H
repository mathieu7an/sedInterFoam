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

if (favreAveraging)
{
    #include "calcFavreAverage.H"
}
if (runTime.outputTime())
{
    if (writeTau == 1)
    {
        Taus = rhos*(alphas*nuEffs+nuFrs)*dev(fvc::grad(Us));
        Tauw = rhow*(scalar(1.0)-alphas)*nuEfff*dev(fvc::grad(Uw));
        Taus.write();
        Tauw.write();
    }
    if (granularRheology.on())
    {
        ps.write();
        granularRheology.muI().write();
        if (granularRheology.dilatancy())
        {
            alphasPlastic.write();
            delta.write();
        }
    }
    if (kineticTheory.on())
    {
        ps.write();
        kineticTheory.Theta().write();
    }
    if (SUSlocal)
    {
        SUS.write();
    }
    pff.write();
    runTime.write();
}
