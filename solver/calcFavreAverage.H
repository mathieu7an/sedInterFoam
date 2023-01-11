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

alphawUwMean = (alphaw*Uw + avCount*alphawUwMean)/(avCount+1);
alphasUsMean = (alphas*Us + avCount*alphasUsMean)/(avCount+1);
alpha_sMean = (alphas + avCount*alpha_sMean)/(avCount+1);

pMeanF = (p_rgh + avCount*pMeanF)/(avCount+1);

UwMeanF = alphawUwMean/(1-alpha_sMean);
UsMeanF = alphasUsMean/(alpha_sMean+1e-20);

UwPrimeF = Uw - UwMeanF;
UsPrimeF = Us - UsMeanF;

pPrimeF = p_rgh - pMeanF;

volScalarField UwPrimeFx = UwPrimeF.component(vector::X);
volScalarField UwPrimeFy = UwPrimeF.component(vector::Y);

dimensionedScalar zeroo
(
"zeroo",
dimensionSet(0, 2, -2, 0, 0, 0, 0),
scalar(0.0)
);

RS2 = min(UwPrimeFx*UwPrimeFy*neg(UwPrimeFx), zeroo);
RS4 = min(UwPrimeFx*UwPrimeFy*pos(UwPrimeFx), zeroo);

RS2Mean = (RS2 + avCount*RS2Mean)/(avCount+1);
RS4Mean = (RS4 + avCount*RS4Mean)/(avCount+1);

alphawSqrUwMean = (alphaw*Uw*Uw + avCount*alphawSqrUwMean)/(avCount+1);
alphasSqrUsMean = (alphas*Us*Us + avCount*alphasSqrUsMean)/(avCount+1);

UwPrime2MeanF = (alphawSqrUwMean/(1-alpha_sMean))-UwMeanF*UwMeanF;
UsPrime2MeanF = alphasSqrUsMean/(alpha_sMean+1e-20)-UsMeanF*UsMeanF;

if (runTime.outputTime())
{
    alpha_sMean.write();
    RS2Mean.write();
    RS4Mean.write();
    UwMeanF.write();
    UsMeanF.write();
    UwPrime2MeanF.write();
    UsPrime2MeanF.write();
}

if (TKEBudgetLES)
{
    Info<<"Compute TKE budget for LES"<<endl;
    TKEMean = 0.5 * (
                       UwPrime2MeanF.component(0)
                     + UwPrime2MeanF.component(3)
                     + UwPrime2MeanF.component(5)
                     );
    TKEMeanProd = -UwPrime2MeanF && fvc::grad(UwMeanF);
    BB = - alphaw*(2.0*symm(fvc::grad(Uw))
    -((2.0/3.0)*I)*tr(fvc::grad(Uw)().T()));

    turbDiffMean = (turbDiffMean - avCount*0.5*(UwPrimeF*magSqr(UwPrimeF)))/
                    (avCount+1);
    turbDiffusionMean = fvc::div(turbDiffMean);
    //Pressure Diffusion Term
    pressDiffusionMean = (pressDiffusionMean - avCount*1./phasew->rho()*
                    (UwPrimeF & fvc::grad(pPrimeF)))/
                    (avCount+1);
    SGSDiffMean = (SGSDiffMean -avCount * (UwPrimeF & (turbulencef->nut()*
                    dev(twoSymm(fvc::grad(Uw))))))/(avCount+1);
    SGSDiffusionMean = fvc::div(SGSDiffMean);
    //Tensor of Strain Rate of Resolved Fluctuations
    SGSstrainTensor = symm(fvc::grad(UwPrimeF));
    //Viscous Dissipation of Resolved Fluctuations
    viscDissMean = (viscDissMean - avCount*2*phasew->nu()*
                   (SGSstrainTensor && SGSstrainTensor))/(avCount+1);
    SGSDissMean = (SGSDissMean + avCount*(turbulencef->nut()*BB &&
                   SGSstrainTensor))/(avCount+1);

    dragContwMean = ((-UwPrimeF & (alphas*alphaw*drags->K(mag(Uw-Us))*
                     (Uw-Us))/rhow)+ avCount*dragContwMean)/(avCount+1);
//    dragContbMean = ((-UwPrimeF & (alpha*alphaw*(Uw-Us))/rhow)
//                     + avCount*dragContbMean)/(avCount+1);

    if (runTime.outputTime())
    {
        TKEMeanProd.write();
        turbDiffusionMean.write();
        pressDiffusionMean.write();
        SGSDiffusionMean.write();
        viscDissMean.write();
        SGSDissMean.write();
        dragContwMean.write();
    }
}

avCount++;
