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

#include "energyJump.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace compressiblePhaseChangeMixtures
{
    defineTypeNameAndDebug(energyJump, 0);
    addToRunTimeSelectionTable
    (
        compressiblePhaseChangeMixture,
        energyJump,
        components
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::compressiblePhaseChangeMixtures::energyJump::energyJump
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
Foam::compressiblePhaseChangeMixtures::energyJump::mDot() const
{
    const volScalarField& T =
        mesh_.lookupObject<volScalarField>("T").oldTime();

    dimensionedScalar L = hf_;
    
    volScalarField limitedAlpha1
    (
        min(max(mixture_.alpha1(), scalar(0)), scalar(1))
    );

    volScalarField limitedAlpha2
    (
        min(max(mixture_.alpha2(), scalar(0)), scalar(1))
    );

    tmp<volScalarField> tkappa = mixture_.kappa();
    volScalarField kappa = tkappa.ref();

    volVectorField gradT(fvc::grad(T));
    volVectorField gradAlpha1(fvc::grad(limitedAlpha1));
    volVectorField gradAlpha2(fvc::grad(limitedAlpha2));

    const dimensionedScalar m0(dimDensity/dimTime, Zero);

    return Pair<tmp<volScalarField>>
    (
        max(coeffE_*kappa*(gradT&gradAlpha1)*limitedAlpha1/L,m0),
        max(coeffC_*kappa*(gradT&gradAlpha2)*limitedAlpha2/L,m0)
    );
}

void Foam::compressiblePhaseChangeMixtures::energyJump::correct()
{
}


bool Foam::compressiblePhaseChangeMixtures::energyJump::read()
{
    if (compressiblePhaseChangeMixture::read())
    {
        return true;
    }

    return false;
}


// ************************************************************************* //
