/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 OpenFOAM Foundation
    Copyright (C) 2016-2018 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "wallHeatFluxTwoPhase.H"
#include "turbulentFluidThermoModel.H"
#include "solidThermo.H"
#include "surfaceInterpolate.H"
#include "fvcSnGrad.H"
#include "wallPolyPatch.H"
#include "turbulentFluidThermoModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(wallHeatFluxTwoPhase, 0);
    addToRunTimeSelectionTable(functionObject, wallHeatFluxTwoPhase, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::wallHeatFluxTwoPhase::writeFileHeader(Ostream& os) const
{
    // Add headers to output data
    writeHeader(os, "Wall heat-flux");
    writeCommented(os, "Time");
    writeTabbed(os, "patch");
    writeTabbed(os, "max");
    writeTabbed(os, "integral");
    writeTabbed(os, "liquid");
    writeTabbed(os, "vapor");
    writeTabbed(os, "total HTC");
    writeTabbed(os, "liquid HTC");
    writeTabbed(os, "vapor HTC");

    os  << endl;
}


void Foam::functionObjects::wallHeatFluxTwoPhase::calcHeatFlux
(
    const volScalarField& alpha,
    const volScalarField& kappa,
    volScalarField& wallHeatFluxTwoPhase
)
{
    volScalarField::Boundary& wallHeatFluxTwoPhaseBf = wallHeatFluxTwoPhase.boundaryFieldRef();

    const volScalarField T = mesh_.lookupObject<volScalarField>("T");
    // - Temperature gradient at boundary.
    volScalarField gradT = mag(fvc::grad(T));
    
    const volScalarField::Boundary& TBf = T.boundaryField();
    const volScalarField::Boundary& gradTBf = gradT.boundaryField();
    // - Thermal conductivity boudnary field.
    const volScalarField::Boundary& kappaBf = kappa.boundaryField();
    const volScalarField alphal = mesh_.lookupObject<volScalarField>("alpha.liquid");

    const scalar TAveLiquid = gSum((T*alphal*mesh_.V())())  //dereference from tmp object.
                            / gSum((alphal*mesh_.V())());

    dimensionedScalar Tave("Tave", dimTemperature, TAveLiquid);
    volScalarField Tref
    (
        IOobject
        (
            "Tref",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(Tave)
    );

    const volScalarField::Boundary& TrefBf = Tref.boundaryField();

    // - Update flux into liquid.
    forAll(wallHeatFluxTwoPhaseBf, patchi)
    {
        wallHeatFluxTwoPhaseBf[patchi] = gradTBf[patchi]*kappaBf[patchi];
    }
    // - Update htc
    volScalarField::Boundary& htcBf = htc_.boundaryFieldRef();

    forAll(htcBf, patchi)
    {
        htcBf[patchi] = wallHeatFluxTwoPhaseBf[patchi] / (TBf[patchi] - TrefBf[patchi]);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::wallHeatFluxTwoPhase::wallHeatFluxTwoPhase
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(obr_, name, typeName, dict),
    patchSet_(),
    qrName_("qr"),
    htc_
    (
        volScalarField
        (
            IOobject
            (
                "liquidHtc",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimPower/dimArea/dimTemperature, Zero)
        )
    )
{
    volScalarField* wallHeatFluxTwoPhasePtr
    (
        new volScalarField
        (
            IOobject
            (
                type(),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimMass/pow3(dimTime), Zero)
        )
    );

    mesh_.objectRegistry::store(wallHeatFluxTwoPhasePtr);

    read(dict);

    writeFileHeader(file());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::wallHeatFluxTwoPhase::~wallHeatFluxTwoPhase()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::wallHeatFluxTwoPhase::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    patchSet_ =
        mesh_.boundaryMesh().patchSet
        (
            wordReList(dict.lookupOrDefault("patches", wordReList()))
        );

    dict.readIfPresent("qr", qrName_);

    Info<< type() << " " << name() << ":" << nl;

    if (patchSet_.empty())
    {
        forAll(pbm, patchi)
        {
            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                patchSet_.insert(patchi);
            }
        }

        Info<< "    processing all wall patches" << nl << endl;
    }
    else
    {
        Info<< "    processing wall patches: " << nl;
        labelHashSet filteredPatchSet;
        for (const label patchi : patchSet_)
        {
            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                filteredPatchSet.insert(patchi);
                Info<< "        " << pbm[patchi].name() << endl;
            }
            else
            {
                WarningInFunction
                    << "Requested wall heat-flux on non-wall boundary "
                    << "type patch: " << pbm[patchi].name() << endl;
            }
        }

        Info<< endl;

        patchSet_ = filteredPatchSet;
    }

    return true;
}


bool Foam::functionObjects::wallHeatFluxTwoPhase::execute()
{
    volScalarField& wallHeatFluxTwoPhase = lookupObjectRef<volScalarField>(type());

    if
    (
        foundObject<compressible::turbulenceModel>
        (
            turbulenceModel::propertiesName
        )
    )
    {
        const compressible::turbulenceModel& turbModel =
            lookupObject<compressible::turbulenceModel>
            (
                turbulenceModel::propertiesName
            );

        calcHeatFlux
        (
            turbModel.alphaEff()(),
            turbModel.kappa(),
            wallHeatFluxTwoPhase
        );
    }
    else if (foundObject<fluidThermo>(fluidThermo::dictName))
    {
        const fluidThermo& thermo =
            lookupObject<fluidThermo>(fluidThermo::dictName);

        calcHeatFlux
        (
            thermo.alpha(),
            thermo.kappa(),
            wallHeatFluxTwoPhase
        );
    }
    else if (foundObject<solidThermo>(solidThermo::dictName))
    {
        const solidThermo& thermo =
            lookupObject<solidThermo>(solidThermo::dictName);

        calcHeatFlux(thermo.alpha(), thermo.kappa(), wallHeatFluxTwoPhase);
    }
    else
    {
        FatalErrorInFunction
            << "Unable to find compressible turbulence model in the "
            << "database" << exit(FatalError);
    }

    return true;
}


bool Foam::functionObjects::wallHeatFluxTwoPhase::write()
{
    const volScalarField& wallHeatFluxTwoPhase = lookupObject<volScalarField>(type());

    Log << type() << " " << name() << " write:" << nl
        << "    writing field " << wallHeatFluxTwoPhase.name() << endl;

    // wallHeatFluxTwoPhase.write();
    const volScalarField alphal = mesh_.lookupObject<volScalarField>("alpha.liquid");

    const fvPatchList& patches = mesh_.boundary();

    const surfaceScalarField::Boundary& magSf =
        mesh_.magSf().boundaryField();

    for (const label patchi : patchSet_)
    {
        const fvPatch& pp = patches[patchi];

        const scalarField& hfp = wallHeatFluxTwoPhase.boundaryField()[patchi];
        const scalarField totalHf = magSf[patchi]*hfp;
        const scalar minHfp = gMin(hfp);
        const scalar maxHfp = gMax(hfp);
        const scalar integralHfp = gSum(totalHf);
        // liquid and vapor phase flux.
        const scalarField& alphalp = alphal.boundaryField()[patchi];
        const scalar liquidHfp = gSum((totalHf*alphalp)());
        const scalar vaporHfp = gSum((totalHf*(1-alphalp))());

        const scalarField& htcp= htc_.boundaryField()[patchi];
        const scalarField totalhtcp = magSf[patchi]*htcp;

        const scalar htcpAve = gSum(totalhtcp)/gSum(magSf[patchi]);
        const scalar liquidHtcpAve = gSum((totalhtcp*alphalp)())/gSum((magSf[patchi]*alphalp)());
        const scalar vaporHtcpAve = gSum((totalhtcp*(1-alphalp))())/gSum((magSf[patchi]*(1-alphalp))());
        if (Pstream::master())
        {
            writeCurrentTime(file());

            file()
                << token::TAB << pp.name()
                << token::TAB << maxHfp
                << token::TAB << integralHfp                
                << token::TAB << liquidHfp
                << token::TAB << vaporHfp
                << token::TAB << htcpAve
                << token::TAB << liquidHtcpAve
                << token::TAB << vaporHtcpAve
                << token::TAB << "debug"
                << token::TAB << min(htcp)
                << token::TAB << min(totalhtcp)
                << token::TAB << min(magSf[patchi])
                << token::TAB << min(alphalp)
                << endl;
        }

        Log << "    min/max/integ(" << pp.name() << ") = "
            << minHfp << ", " << maxHfp << ", " << integralHfp << endl;
    }

    return true;
}


// ************************************************************************* //
