/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
    Copyright (C) YEAR AUTHOR, AFFILIATION
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

#include "zeroMassFluxFvPatchFieldTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"
#include "PatchFunction1.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makeRemovablePatchTypeField
(
    fvPatchScalarField,
    zeroMassFluxFvPatchScalarField
);

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::
zeroMassFluxFvPatchScalarField::
zeroMassFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    parent_bctype(p, iF)
{
    if (false)
    {
        //printMessage("Construct zeroMassFlux : patch/DimensionedField");
    }
}


Foam::
zeroMassFluxFvPatchScalarField::
zeroMassFluxFvPatchScalarField
(
    const zeroMassFluxFvPatchScalarField& rhs,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    parent_bctype(rhs, p, iF, mapper)
{
    if (false)
    {
        //printMessage("Construct zeroMassFlux : patch/DimensionedField/mapper");
    }
}


Foam::
zeroMassFluxFvPatchScalarField::
zeroMassFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    parent_bctype(p, iF, dict)
{
    if (false)
    {
        //printMessage("Construct zeroMassFlux : patch/dictionary");
    }
}


Foam::
zeroMassFluxFvPatchScalarField::
zeroMassFluxFvPatchScalarField
(
    const zeroMassFluxFvPatchScalarField& rhs
)
:
    parent_bctype(rhs),
    dictionaryContent(rhs)
{
    if (false)
    {
        //printMessage("Copy construct zeroMassFlux");
    }
}


Foam::
zeroMassFluxFvPatchScalarField::
zeroMassFluxFvPatchScalarField
(
    const zeroMassFluxFvPatchScalarField& rhs,
    const DimensionedField<scalar, volMesh>& iF
)
:
    parent_bctype(rhs, iF)
{
    if (false)
    {
        //printMessage("Construct zeroMassFlux : copy/DimensionedField");
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::
zeroMassFluxFvPatchScalarField::
~zeroMassFluxFvPatchScalarField()
{
    if (false)
    {
        //printMessage("Destroy zeroMassFlux");
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void
Foam::
zeroMassFluxFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (false)
    {
        //printMessage("updateCoeffs zeroMassFlux");
    }

//{{{ begin code
    //#line 31 "/home/thiago/OpenFOAM/v2206/OpenFOAM-v2206/tutorials/BEI/dboukFoam/buoyant_suspension/fcteCase/nouveauC/0/C.boundaryField.movingWall"
// dictionnaires
            const dictionary& transportProperties = db().lookupObject<IOdictionary> ("transportProperties");
            dimensionedScalar DC("DC",dimViscosity,transportProperties);
            dimensionedScalar Cmax("Cmax",transportProperties);
            dimensionedVector vst("vst",transportProperties);

            // Grandeurs "frontieres"
            scalarField& Cf = *this; // concentration

            // Grandeurs "centres"
            const tmp<vectorField>& tvn = patch().nf();
            const vectorField& vn = tvn();

            const tmp<scalarField>& tdelta = this->patch().deltaCoeffs();
            const scalarField& delta = tdelta();
            const tmp<scalarField>& Cc = patchInternalField(); // concentration
            
            
            forAll(Cf, faceID)
            {
                Cf[faceID] = Cc.ref()[faceID]*(1.+(vst.value()&vn[faceID])/(DC.value()*delta[faceID]));
            }

    this->parent_bctype::updateCoeffs();
}


// ************************************************************************* //

