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
// Compute two-phase specific terms for RAS turbulence models:
// drag dissipation and stratification terms.

// Exponent of the turbulence damping term for k,epsilon and omega equations
tmfexp = -B*rhos*6.0/(max(alphaf*K, Tpsmall));

// just need the horizontal ones;
rtemp = (Uf-Us) & gravity;
Urtemp = Uf-Us;
if (Foam::mag(gravity).value() != 0)
{
    Urtemp = rtemp*gravity/magSqr(gravity);
}
Urh = (Uf-Us) - Urtemp;

if (not SUSaniso)
{
// ESD1 for horizontal, ESD4 for vertical density stratification
// for 1D ESD1 ==0, obviously
ESD1 = (K/rhof)*(fvc::grad(alphas)) & Urh;

// stratification
ESD4 = -(SUS/alphaf)*(rhos/rhof-1.0)*gravity & fvc::grad(alphas);
ESD4 = ESD4 + KE1*ESD1;

// the following term is implemented in epsilon equation.
ESD5 = -(SUS/alphaf)*(rhos/rhof-1.0)*gravity & fvc::grad(alphas);
}
else
{
// ESD1 for horizontal, ESD4 for vertical density stratification
// for 1D ESD1 ==0, obviously
ESD1 = SUS*(K/rhof)*(SUS_I*iso-SUS_A*aniso)&(fvc::grad(alphas))
          & Urh;

// stratification
ESD4 = - SUS/alphaf*(rhos/rhof-1.0)*gravity
            & fvc::grad(alphas);
ESD4 = ESD4 + KE1*ESD1;

// the following term is implemented in epsilon equation.
ESD5 = - SUS/alphaf*(rhos/rhof-1.0)*gravity
            & fvc::grad(alphas);
}
// if ESD5 > 0, then it should be implemented, else, no
ESD5 = max(ESD5, scalar(0.0)*ESD5);
ESD5 = ESD5 + KE1*ESD1;

//         k was included in the equation implicitly
ESD3 = -2.0*K*alphas/rhof;

ESD = KE3*ESD3;   //need to check here!
