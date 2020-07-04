
solidThermo & thermo = sthermoPtr();

// const volScalarField& kappa = thermo.kappa();
// const volScalarField& rho = thermo.rho();
// const volScalarField& Cp = thermo.Cp();
tmp<volScalarField> trho = thermo.rho();
const volScalarField& rho = trho();


// volScalarField& T = thermo.T();
volScalarField& h = thermo.he();


// volScalarField DT = kappa/(rho*Cp);