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
    pcSu_
    (
        volScalarField::Internal
        (
            IOobject
            (
                "pcSu",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(dimless/dimTime, Zero)
        )
    ),
    pcSp_
    (
        volScalarField::Internal
        (
            IOobject
            (
                "pcSp",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(dimless/dimTime, Zero)
        )
    ),
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
    T0_
    (
        "T0",dimTemperature,subDict("saturationProperty")
    ),
    TbyP_
    (
        "TbyP",dimTemperature/dimPressure,subDict("saturationProperty")
    )
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

    Pair<tmp<volScalarField>> mDot = this->mDot();
    const volScalarField& mDot1 = mDot[0]();
    const volScalarField& mDot2 = mDot[1]();

    const volScalarField dmdtNet(mDot1 - mDot2);
    Info<<"dmdtNet Min/Max"<< min(dmdtNet).value() << ", " << max(dmdtNet).value() << endl;

    volScalarField heat_source(dmdtNet * hf_ / mixture_.Cv());
    Info<<"T source Min/Max"<< min(heat_source).value() << ", " << max(heat_source).value() << endl;

    eqn -=
    (
        dmdtNet * hf_ / mixture_.Cv()
    );

    return tEqnPtr;
}

Foam::tmp<Foam::volScalarField> Foam::compressiblePhaseChangeMixture::dmdtNet()
{
    Pair<tmp<volScalarField>> mDot = this->mDot();
    const volScalarField& mDot1 = mDot[0]();
    const volScalarField& mDot2 = mDot[1]();
    
    dmdtNet_ = mDot2 -mDot1;

    return dmdtNet_;
}


Foam::dimensionedScalar Foam::compressiblePhaseChangeMixture::updatedTsat() const
{
    const volScalarField& p = mesh_.lookupObject<volScalarField>("p_rgh");

    dimensionedScalar pmax = max(p);

    const dimensionedScalar Tsatp(T0_ + pmax*TbyP_);

    return Tsatp;
}

Foam::tmp<Foam::volScalarField> 
Foam::compressiblePhaseChangeMixture::coeffs() const
{
    return (1.0/mixture_.thermo1().rho()-1.0/mixture_.thermo2().rho());
}

Foam::volScalarField::Internal & Foam::compressiblePhaseChangeMixture::pcSu()
{
    return pcSu_;
}

Foam::volScalarField::Internal & Foam::compressiblePhaseChangeMixture::pcSp()
{
    return pcSp_;
}


Foam::tmp<Foam::volScalarField> 
Foam::compressiblePhaseChangeMixture::rho1()
{
    return mixture_.thermo1().rho();
}

Foam::tmp<Foam::volScalarField> 
Foam::compressiblePhaseChangeMixture::rho2()
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
