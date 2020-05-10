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
// NOTE: dmdtNet is distributed in terms =
//  Source for phase 1 =
//      dmdtNet/rho1
//    - alpha1*dmdtNet(1/rho1 - 1/rho2)
//  second term is added to Su or Sp regarding the 
//  sing of dmdtNet

tmp<volScalarField> tdmdtNet = fluid.dmdtNet();
const volScalarField& dmdt12 = tdmdtNet();

tmp<volScalarField> tcoeffs(fluid.coeffs());
const volScalarField& coeffs = tcoeffs();

forAll(dmdt12, celli)
{
    scalar cdmdt = dmdt12[celli];
    scalar ccoeffs = coeffs[celli];

    scalar alpha1Limited = max(min(alpha1[celli],1.0),0.0);

    Su[celli] += 1.0/rho1[celli]*cdmdt;

    if (cdmdt > 0)
    {
        Su[celli] -= cdmdt * ccoeffs * alpha1Limited;
    }
    else if (cdmdt < 0)
    {
        Sp[celli] -= cdmdt * ccoeffs;
    }
}

// // - All terms are treated explicitly
// Su +=  dmdt12/rho1 - alpha1 * dmdt12* fluid.coeffs();



volScalarField::Internal divU
(
    mesh.moving()
  ? fvc::div(phiCN() + mesh.phi())
  : fvc::div(phiCN())
);
