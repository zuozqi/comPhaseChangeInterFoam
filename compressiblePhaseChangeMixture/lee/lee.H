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
    Foam::phaseChangeTwoPhaseMixtures::lee

Description
    lee condensation/saturation model.


SourceFiles
    lee.C

\*--------------------------------------------------------------------*/

#ifndef lee_H
#define lee_H

#include "compressiblePhaseChangeMixture.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressiblePhaseChangeMixtures
{

/*--------------------------------------------------------------------*\
                             Class lee
\*--------------------------------------------------------------------*/

class lee
:
    public compressiblePhaseChangeMixture
{
    // Private data

        //- Condensation coefficient [1/s/K]
        dimensionedScalar coeffC_;

        //- Evaporation coefficient [1/s/K]
        dimensionedScalar coeffE_;


public:

    //- Runtime type information
    TypeName("lee");


    // Constructors

        //- Construct from components
        lee
        (
            const twoPhaseMixtureThermo& mixture,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~lee() = default;


    // Member Functions
        //- Return the mass condensation and vaporisation rates as coefficients
        virtual Pair<tmp<volScalarField>> mDot();
        //- Correct the lee phaseChange model
        virtual void correct();

        //- Read the transportProperties dictionary and update
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressiblePhaseChangeMixtures
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
