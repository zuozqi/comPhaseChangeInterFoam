/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2019 OpenCFD Ltd.
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

#include "compressiblePhaseChangeMixture.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(compressiblePhaseChangeMixture, 0);
    defineRunTimeSelectionTable
    (
        compressiblePhaseChangeMixture,
        components
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::compressiblePhaseChangeMixture::
compressiblePhaseChangeMixture
(
    const twoPhaseMixtureThermo& mixture,
    const fvMesh& mesh
)
:
    IOdictionary
    (
        IOobject
        (
            "phaseChangeProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    mixture_(mixture),
    mesh_(mesh),
    dmdtNet_
    (
        volScalarField
        (
            IOobject
            (
                "dmdtNet",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar(dimDensity/dimTime, Zero)
        )
    ),
    hf_
    (
        "hf",dimEnergy/dimMass,subDict("saturationProperty")
    ),
    C_(subDict("saturationProperty").lookup("C<8>"))
{

}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::fvScalarMatrix> 
Foam::compressiblePhaseChangeMixture::heatTransfer(const volScalarField & T)
{
    Foam::tmp<Foam::fvScalarMatrix> tEqnPtr
    (
        new fvScalarMatrix(T, dimMass*dimTemperature/dimTime)
    );

    fvScalarMatrix& eqn = tEqnPtr.ref();


    tmp<volScalarField> tdmdt12 = dmdtNet();
    const volScalarField dmdt12 = tdmdt12();
    Info<<"dmdtNet Min/Max = "<< min(dmdt12).value() << ", " << max(dmdt12).value() << endl;

    eqn +=
    (
        dmdt12 * hf_ / mixture_.Cv()
    );

    return tEqnPtr;
}

Foam::tmp<Foam::volScalarField> Foam::compressiblePhaseChangeMixture::dmdtNet()
{
    Pair<tmp<volScalarField>> mDot = this->mDot();
    const volScalarField& mDot1 = mDot[0]();
    const volScalarField& mDot2 = mDot[1]();

    dmdtNet_ = mDot2 - mDot1;

    return dmdtNet_;
}


Foam::tmp<Foam::volScalarField>
Foam::compressiblePhaseChangeMixture::updatedTsat() const
{
    const volScalarField& p = mesh_.lookupObject<volScalarField>("p_rgh");
    tmp<volScalarField> tTsat
    (
        volScalarField::New
        (
            "Tsat",
            p.mesh(),
            dimensionedScalar(dimTemperature)
        )
    );

    volScalarField& Tsat = tTsat.ref();

    forAll(Tsat, celli)
    {
        Tsat[celli] = C_.value(p[celli]);
    }

    volScalarField::Boundary& TsatBf = Tsat.boundaryFieldRef();

    forAll(Tsat.boundaryField(), patchi)
    {
        scalarField& Tsatp = TsatBf[patchi];
        const scalarField& pp = p.boundaryField()[patchi];

        forAll(Tsatp, facei)
        {
            Tsatp[facei] = C_.value(pp[facei]);
        }
    }

    return tTsat;
}

Foam::tmp<Foam::volScalarField> 
Foam::compressiblePhaseChangeMixture::coeffs() const
{
    return (1.0/mixture_.thermo1().rho()-1.0/mixture_.thermo2().rho());
}

Foam::tmp<Foam::volScalarField> 
Foam::compressiblePhaseChangeMixture::rho1() const
{
    return mixture_.thermo1().rho();
}

Foam::tmp<Foam::volScalarField> 
Foam::compressiblePhaseChangeMixture::rho2() const
{
    return mixture_.thermo2().rho();
}

bool Foam::compressiblePhaseChangeMixture::read()
{
    if (regIOobject::read())
    {
        return true;
    }

    return false;
}


// ************************************************************************* //
