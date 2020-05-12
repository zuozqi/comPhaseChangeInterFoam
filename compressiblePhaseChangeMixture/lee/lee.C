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

#include "lee.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace compressiblePhaseChangeMixtures
{
    defineTypeNameAndDebug(lee, 0);
    addToRunTimeSelectionTable
    (
        compressiblePhaseChangeMixture,
        lee,
        components
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::compressiblePhaseChangeMixtures::lee::lee
(
    const twoPhaseMixtureThermo& mixture,
    const fvMesh& mesh
)
:
    compressiblePhaseChangeMixture(mixture, mesh),
    coeffC_
    (
        "coeffC", dimless/dimTime/dimTemperature, subDict("Coeffs")
    ),
    coeffE_
    (
        "coeffE", dimless/dimTime/dimTemperature, subDict("Coeffs")
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::compressiblePhaseChangeMixtures::lee::mDot() const
{

    volScalarField limitedAlpha1
    (
        min(max(mixture_.alpha1(), scalar(0)), scalar(1))
    );

    volScalarField limitedAlpha2
    (
        min(max(mixture_.alpha2(), scalar(0)), scalar(1))
    );

    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");

    const volScalarField& rho1 = mixture_.thermo1().rho();
    const volScalarField& rho2 = mixture_.thermo2().rho();

    tmp<volScalarField> tTsat = this->updatedTsat();
    const volScalarField Tsat = tTsat();
    
    const dimensionedScalar T0(dimTemperature, Zero);

    return Pair<tmp<volScalarField>>
    (
        coeffE_*rho1*limitedAlpha1*max(T.oldTime() - Tsat, T0),
        coeffC_*rho2*limitedAlpha2*max(Tsat - T.oldTime(), T0)
    );
}

void Foam::compressiblePhaseChangeMixtures::lee::correct()
{
}


bool Foam::compressiblePhaseChangeMixtures::lee::read()
{
    if (compressiblePhaseChangeMixture::read())
    {
        subDict("Coeffs").readEntry("coeffC", coeffC_);
        subDict("Coeffs").readEntry("coeffE", coeffE_);

        return true;
    }

    return false;
}


// ************************************************************************* //
