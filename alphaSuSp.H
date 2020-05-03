volScalarField::Internal Sp
(
    IOobject
    (
        "Sp",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(dgdt.dimensions(), Zero)
);

volScalarField::Internal Su
(
    IOobject
    (
        "Su",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(dgdt.dimensions(), Zero)
);

forAll(dgdt, celli)
{
    if (dgdt[celli] > 0.0)
    {
        Sp[celli] -= dgdt[celli]/max(1.0 - alpha1[celli], 1e-4);
        Su[celli] += dgdt[celli]/max(1.0 - alpha1[celli], 1e-4);
    }
    else if (dgdt[celli] < 0.0)
    {
        Sp[celli] += dgdt[celli]/max(alpha1[celli], 1e-4);
    }
}

tmp<volScalarField> tdmdtNet = fluid.dmdtNet();
const volScalarField& dmdt12 = tdmdtNet();
        // NOTE: dmdtNet is distributed in terms =
        //  Source for phase 1 =
        //      dmdtNet/rho1
        //    - alpha1*dmdtNet(1/rho1 - 1/rho2)

Su +=  dmdt12/rho1 - alpha1 * dmdt12* fluid.coeffs();
// fluid.updateSuSp();
// volScalarField::Internal& pcSu(fluid.pcSu());
// volScalarField::Internal& pcSp(fluid.pcSp());

// //- update phase change source term.
// Su += pcSu;
// Sp += pcSp;


volScalarField::Internal divU
(
    mesh.moving()
  ? fvc::div(phiCN() + mesh.phi())
  : fvc::div(phiCN())
);
