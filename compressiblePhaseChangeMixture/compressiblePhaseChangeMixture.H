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

Class
    Foam::compressiblePhaseChangeMixture

Description

SourceFiles
    compressiblePhaseChangeMixture.C
    newcompressiblePhaseChangeMixture.C

\*---------------------------------------------------------------------------*/

#ifndef compressiblePhaseChangeMixture_H
#define compressiblePhaseChangeMixture_H

#include "twoPhaseMixtureThermo.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "fvScalarMatrix.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class compressiblePhaseChangeMixture Declaration
\*---------------------------------------------------------------------------*/

class compressiblePhaseChangeMixture
:
    public IOdictionary
{

protected:

    // Protected data

        //- Reference to the twoPhaseMixtureThermo
        const twoPhaseMixtureThermo& mixture_;

        //- Reference to fvMesh
        const fvMesh& mesh_;

        volScalarField::Internal pcSu_;
        //- Sp phase source term
        volScalarField::Internal pcSp_;

        //- Mass flux, kg/m^3/s.
        volScalarField dmdtNet_;

        //- Saturation properties
        dimensionedScalar hf_;
        dimensionedScalar T0_;     // For nitrogen: T0:69.954, TbyP:7.32646e-5;
        dimensionedScalar TbyP_;  //- T:K P:Pa Tsatp = T0 + Pmax*TbyP;

    // Private Member Functions

        //- No copy construct
        compressiblePhaseChangeMixture
        (
            const compressiblePhaseChangeMixture&
        ) = delete;

        //- No copy assignment
        void operator=(const compressiblePhaseChangeMixture&) = delete;


public:

    //- Runtime type information
    TypeName("compressiblePhaseChangeMixture");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            compressiblePhaseChangeMixture,
            components,
            (
                const twoPhaseMixtureThermo& mixture,
                const fvMesh& mesh
            ),
            (mixture, mesh)
        );


    // Selectors

        //- Return a reference to the selected phaseChange model
        static autoPtr<compressiblePhaseChangeMixture> New
        (
            const twoPhaseMixtureThermo& mixture,
            const fvMesh& mesh
        );


    // Constructors

        //- Construct from components
        compressiblePhaseChangeMixture
        (
            const twoPhaseMixtureThermo& mixture,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~compressiblePhaseChangeMixture() = default;


    // Member Functions

        //- The mass flux rate of evaporation/condensation, in unit of kg/m^3/s
        virtual Pair<tmp<volScalarField>> mDot() const = 0;

        //- Net mass transfer from phase1(liquid) to phase2(vapor)
        tmp<volScalarField> dmdtNet();

        //Return the coefficients in PEqn and alphaEqn, (1/rho1 -1/rho2).
        virtual tmp<volScalarField> coeffs() const;

        // //- update the Su and Sp according to the phase change rate.
        // void updateSuSp();
        
        //- updated saturation temperature according to polynomial fitting

        dimensionedScalar updatedTsat() const;

        // Calculate the energy change due to phase chagne.
        virtual tmp<fvScalarMatrix> heatTransfer(const volScalarField& T);


        volScalarField::Internal & pcSu();
        volScalarField::Internal & pcSp();

        //- Correct the phaseChange model
        virtual void correct() = 0;

        virtual tmp<volScalarField> rho1();
        virtual tmp<volScalarField> rho2();
        //- Read the transportProperties dictionary and update
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
