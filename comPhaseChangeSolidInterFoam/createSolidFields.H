// Info<< "Reading solid field sT\n" << endl;

// volScalarField T
// (
//     IOobject
//     (
//         "T",
//         runTime.timeName(),
//         solidMesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE
//     ),
//     solidMesh,
//     dimensionedScalar(dimTemperature, Zero)
// );

Info<< "Creating solid Thermo"<<endl;

// PtrList<solidThermo> thermos(1);
// thermos.set(0, solidThermo::New(solidMesh));
// solidThermo& sthermo = thermos[0];

autoPtr<solidThermo> sthermoPtr = solidThermo::New(solidMesh);
// solidThermo & sthermo = sthermoPtr();

// const volScalarField& skappa = sthermo.kappa();
// const volScalarField& srho = sthermo.rho();
// const volScalarField& sCp = sthermo.Cp();
// volScalarField& sT = sthermo.T();

// volScalarField DT = skappa/(srho*sCp);



// Info<< "Reading solid diffusivity DT\n" << endl;

// autoPtr<volScalarField> DT
// (
//     new volScalarField
//     (
//     IOobject
//     (
//         "DT",
//         runTime.timeName(),
//         solidMesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE
//     ),
//     solidMesh,
//     dimensionedScalar(dimViscosity, Zero)
//     )
// );
// volScalarField DT
// (
//     IOobject
//     (
//         "DT",
//         runTime.timeName(),
//         solidMesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE
//     ),
//     solidMesh,
//     dimensionedScalar(dimViscosity, Zero)
// );

// // Also can be provided in transportProperties dict located in constant/
// if (!DT.headerOk())
// {
//     IOdictionary transportProperties
//     (
//         IOobject
//         (
//             "transportProperties",
//             runTime.constant(),
//             mesh,
//             IOobject::MUST_READ_IF_MODIFIED,
//             IOobject::NO_WRITE
//         )
//     );
//     DT = dimensionedScalar("DT", dimViscosity, transportProperties);
// }

// #include "createFvOptions.H"
