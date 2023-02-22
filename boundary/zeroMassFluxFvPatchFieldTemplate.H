/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

Description
    Template for use with dynamic code generation of a
    zeroMassFlux fvPatchField.

    - without state

SourceFiles
    zeroMassFluxFvPatchFieldTemplate.C

\*---------------------------------------------------------------------------*/

#ifndef zeroMassFluxFvPatchTemplateScalarField_H
#define zeroMassFluxFvPatchTemplateScalarField_H

#include "fixedValueFvPatchFields.H"
#include "dictionaryContent.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         A templated FixedValueFvPatch
\*---------------------------------------------------------------------------*/

class zeroMassFluxFvPatchScalarField
:
    public fixedValueFvPatchField<scalar>,
    public dictionaryContent
{
    //- The parent boundary condition type
    typedef fixedValueFvPatchField<scalar> parent_bctype;


    // Private Member Functions

        //- Report a message with the SHA1sum
        /*inline static void printMessage(const char* message)
        {
            Info<< message << " sha1: " << SHA1sum << '\n';
        }*/

public:

    //- SHA1 representation of the code content
    //static constexpr const char* const SHA1sum = "a35f374cbca51d51099bcf7fc0521d3c8dccf7d7";

    //- Runtime type information
    TypeName("zeroMassFlux");


    // Constructors

        //- Construct from patch and internal field
        zeroMassFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        zeroMassFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping a copy onto a new patch
        zeroMassFluxFvPatchScalarField
        (
            const zeroMassFluxFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        zeroMassFluxFvPatchScalarField
        (
            const zeroMassFluxFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new zeroMassFluxFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        zeroMassFluxFvPatchScalarField
        (
            const zeroMassFluxFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new zeroMassFluxFvPatchScalarField(*this, iF)
            );
        }


    //- Destructor
    virtual ~zeroMassFluxFvPatchScalarField();


    // Member Functions

        //- Code context as a dictionary
        const dictionary& codeContext() const noexcept
        {
            return dictionaryContent::dict();
        }

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

