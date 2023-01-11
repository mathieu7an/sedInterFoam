
// Cell gradient of alphaw
volVectorField gradAlphaw(fvc::grad(alphaw));;

// Interpolated face-gradient of alphaw
surfaceVectorField gradAlphawf(fvc::interpolate(gradAlphaw));

// Face unit interface normal
surfaceVectorField  nHatfv(gradAlphawf/(mag(gradAlphawf) + deltaN));

const volScalarField::Boundary& alphawbound = alphaw.boundaryField();
//const volScalarField::Boundary& alphasbound = alphas.boundaryField();
const volScalarField::Boundary& alphagbound = alphag.boundaryField();
const volVectorField::Boundary& Ubound = U.boundaryField();

const fvBoundaryMesh& boundary = mesh.boundary();

forAll(boundary, patchi)
{
    if (isA<alphaContactAngleTwoPhaseFvPatchScalarField>(alphawbound[patchi]))
    {
        const alphaContactAngleTwoPhaseFvPatchScalarField& acap =
            refCast<const alphaContactAngleTwoPhaseFvPatchScalarField>
            (alphagbound[patchi]);

        fvsPatchVectorField& nHatp = nHatfv.boundaryFieldRef()[patchi];

        scalarField theta1
        (
            degToRad() * (180 - acap.theta(Ubound[patchi], nHatp))
        );
        
        vectorField nf
        (
            boundary[patchi].nf()
        );

        // Reset nHatPatch to correspond to the contact angle

        scalarField a12(nHatp & nf);

        scalarField b1(cos(theta1));

        scalarField b2(nHatp.size());

        forAll(b2, facei)
        {
            b2[facei] = Foam::cos(Foam::acos(a12[facei]) - theta1[facei]);
        }

        scalarField det(1.0 - a12*a12);

        scalarField a((b1 - a12*b2)/det);
        scalarField b((b2 - a12*b1)/det);

        nHatp = a*nf + b*nHatp;
        nHatp /= (mag(nHatp) + deltaN.value());
    }
}

nHatf = nHatfv & mesh.Sf();

Kg = -fvc::div(nHatf);

//volScalarField limitedAlphas(max(alphas, scalar(0)));
volScalarField limitedAlphag(max(alphag, scalar(0)));

sigma = (limitedAlphag*sigmagf)/(limitedAlphag + SMALL);
