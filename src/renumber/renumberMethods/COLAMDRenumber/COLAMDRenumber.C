/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020-2022 OpenCFD Ltd.
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

#include "COLAMDRenumber.H"
#include "addToRunTimeSelectionTable.H"
#include "bandCompression.H"
#include "decompositionMethod.H"

#include <iostream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(COLAMDRenumber, 0);

    addToRunTimeSelectionTable
    (
        renumberMethod,
        COLAMDRenumber,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::COLAMDRenumber::COLAMDRenumber(const dictionary& dict)
:
    renumberMethod(dict)
{
    std::cout << "hiiiiiiiiii" << std::endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::COLAMDRenumber::renumber
(
    const polyMesh& mesh,
    const pointField& points
) const
{
    CompactListList<label> cellCells;
    decompositionMethod::calcCellCells
    (
        mesh,
        identity(mesh.nCells()),
        mesh.nCells(),
        false,                      // local only
        cellCells
    );

    labelList orderedToOld = meshTools::bandCompression(cellCells);

    return orderedToOld;
}


Foam::labelList Foam::COLAMDRenumber::renumber
(
    const labelList& cellCells,
    const labelList& offsets,
    const pointField& cc
) const
{
    labelList orderedToOld = meshTools::bandCompression(cellCells, offsets);

    return orderedToOld;
}


Foam::labelList Foam::COLAMDRenumber::renumber
(
    const CompactListList<label>& cellCells,
    const pointField& cc
) const
{
    labelList orderedToOld = meshTools::bandCompression(cellCells);

    return orderedToOld;
}


Foam::labelList Foam::COLAMDRenumber::renumber
(
    const labelListList& cellCells,
    const pointField& points
) const
{
    labelList orderedToOld = meshTools::bandCompression(cellCells);

    return orderedToOld;
}


// ************************************************************************* //
