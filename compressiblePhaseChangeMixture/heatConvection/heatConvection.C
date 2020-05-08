/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenCFD Ltd.
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

#include "heatConvection.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace compressiblePhaseChangeMixtures
{
    defineTypeNameAndDebug(heatConvection, 0);
    addToRunTimeSelectionTable
    (
        compressiblePhaseChangeMixture,
        heatConvection,
        components
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::compressiblePhaseChangeMixtures::heatConvection::heatConvection
(
    const twoPhaseMixtureThermo& mixture,
    const fvMesh& mesh
)
:
    compressiblePhaseChangeMixture(mixture, mesh),
    coeffC_
    (
        "coeffC", dimless, subDict("Coeffs")
    ),
    coeffE_
    (
        "coeffE", dimless, subDict("Coeffs")
    )
{

}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::compressiblePhaseChangeMixtures::heatConvection::mDot() const
{
    volScalarField limitedAlpha1
    (
        min(max(mixture_.alpha1(), scalar(0)), scalar(1))
    );

    volScalarField limitedAlpha2
    (
        min(max(mixture_.alpha2(), scalar(0)), scalar(1))
    );

    const volScalarField& T =
        mesh_.lookupObject<volScalarField>("T").oldTime();

    dimensionedScalar L = hf_;

    const dimensionedScalar Tsat = this->updatedTsat();

    const volScalarField areaDensity = mag(fvc::grad(limitedAlpha1));

    const dimensionedScalar Nl
    (
       sum(mesh_.V())/sum(mesh_.V()*areaDensity)
    );

    // <<Heat Transfer>> P292 Eqn 7-31
    dimensionedScalar Ptl
    (
        "Ptl", dimless, subDict("Coeffs")
    );

    dimensionedScalar Ptg
    (
        "Ptg", dimless, subDict("Coeffs")
    );

    volScalarField Pt = Ptl * limitedAlpha1 + Ptg * limitedAlpha2;
    dimensionedScalar g("g",dimVelocity/dimTime,9.81);

    volScalarField beta = scalar(1)/T;

    dimensionedScalar qw
    (
        "qw", dimEnergy/dimLength/dimLength/dimTime, subDict("Coeffs")
    );

    tmp<volScalarField> tkappa = mixture_.kappa();
    volScalarField kappa = tkappa.ref();

    tmp<volScalarField> tnu = mixture_.nu();
    volScalarField nu = tnu.ref();

    volScalarField Grxstar = (g*beta*qw*pow(Nl,4))/(kappa*nu*nu);

    const volScalarField h = 0.60*pow((Grxstar*Pt),0.2)*kappa/Nl;

    const dimensionedScalar T0(dimTemperature, Zero);

    return Pair<tmp<volScalarField>>
    (
        h*max(T-Tsat,T0)*limitedAlpha1/L/Nl,
        h*max(Tsat-T,T0)*limitedAlpha2/L/Nl
    );
}

void Foam::compressiblePhaseChangeMixtures::heatConvection::correct()
{
}


bool Foam::compressiblePhaseChangeMixtures::heatConvection::read()
{
    if (compressiblePhaseChangeMixture::read())
    {
        return true;
    }

    return false;
}


// ************************************************************************* //
