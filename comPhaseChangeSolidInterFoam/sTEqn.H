fvScalarMatrix TEqn
(
        fvm::ddt(rho, h)
        - 
        fvm::laplacian(thermo.alpha(), h, "laplacian(alpha,h)")
        ==
        fvOptions(rho, h)
);
TEqn.relax();

fvOptions.constrain(TEqn);

TEqn.solve();

thermo.correct();

Info<< "Solid Min/max T:" << min(thermo.T()).value() << ' '
    << max(thermo.T()).value() << endl;

fvOptions.correct(h);

