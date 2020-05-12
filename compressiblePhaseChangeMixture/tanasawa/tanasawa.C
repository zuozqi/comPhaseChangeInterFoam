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

#include "constants.H"
#include "tanasawa.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvcDiv.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace compressiblePhaseChangeMixtures
{
    defineTypeNameAndDebug(tanasawa, 0);
    addToRunTimeSelectionTable
    (
        compressiblePhaseChangeMixture,
        tanasawa,
        components
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::compressiblePhaseChangeMixtures::tanasawa::tanasawa
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
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::compressiblePhaseChangeMixtures::tanasawa::mDot() const
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

    tmp<volScalarField> tMv = mixture_.W();
    const volScalarField& Mv = tMv();
    
    dimensionedScalar L = hf_;

    tmp<volScalarField> tTsat = this->updatedTsat();
    const volScalarField Tsat = tTsat();

    tmp<volScalarField> tRhom
    (
        new volScalarField
        (
            IOobject
            (
                "trhom",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimDensity, Zero)
        )
    );

    volScalarField& rhom = tRhom.ref();

    rhom =
        this->rho2()*this->rho1()
        / (this->rho1() - this->rho2());

    const volScalarField tanaConst
    (
        sqrt
        (
            Mv
            /2.0
            /constant::physicoChemical::R
            /constant::mathematical::pi
            /pow3(Tsat)
        )*rhom*L
    );

    const volScalarField areaDensity = mag(fvc::grad(limitedAlpha1));

    const dimensionedScalar Nl
    (
        gSum((areaDensity*mesh_.V())())
        /(
            gSum
            (
                ((areaDensity*limitedAlpha1)*mesh_.V())()
            )
            + dimensionedScalar("SMALL", dimless, VSMALL)
        )
    );

    const dimensionedScalar T0(dimTemperature, Zero);


    return Pair<tmp<volScalarField>>
    (
        2*coeffE_*tanaConst*areaDensity*Nl*limitedAlpha1*(max(T.oldTime() - Tsat, T0))/(2-coeffE_),
        2*coeffC_*tanaConst*areaDensity*Nl*limitedAlpha2*(max(Tsat - T.oldTime(), T0))/(2-coeffC_)
    );
}

void Foam::compressiblePhaseChangeMixtures::tanasawa::correct()
{
}


bool Foam::compressiblePhaseChangeMixtures::tanasawa::read()
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
