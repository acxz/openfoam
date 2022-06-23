/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::metisRenumber

Description
    Renumber using METIS's nested dissection algorithm

SourceFiles
    metisRenumber.C

\*---------------------------------------------------------------------------*/

#include "metisRenumber.H"
#include "addToRunTimeSelectionTable.H"
#include "decompositionMethod.H"
#include "metis.h"
#include "polyMesh.H"
#include "PrecisionAdaptor.H"

#include <iostream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(metisRenumber, 0);

    addToRunTimeSelectionTable
    (
        renumberMethod,
        metisRenumber,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::metisRenumber::metisRenumber(const dictionary& dict)
:
    renumberMethod(dict),
    coeffsDict_(dict.optionalSubDict(typeName+"Coeffs"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::metisRenumber::renumber
(
    const polyMesh& mesh,
    const pointField& points
) const
{

    std::cout << "tryna metis renumber" << std::endl;

    // Make Metis CSR (Compressed Storage Format) storage
    //   adjncy      : contains neighbours (= edges in graph)
    //   xadj(celli) : start of information in adjncy for celli

    CompactListList<label> cellCells;
    decompositionMethod::calcCellCells
    (
        mesh,
        identity(mesh.nCells()),
        mesh.nCells(),
        true,
        cellCells
    );

    const labelList& adjncy = cellCells.m();
    const labelList& xadj = cellCells.offsets();

    idx_t numCells = max(0, (xadj.size()-1));

    // Addressing
    ConstPrecisionAdaptor<idx_t, label, List> xadj_param(xadj);
    ConstPrecisionAdaptor<idx_t, label, List> adjncy_param(adjncy);

    // Avoid potential nullptr issues with zero-sized arrays
    labelList adjncy_dummy, xadj_dummy;
    if (!numCells)
    {
        adjncy_dummy.resize(1, 0);
        adjncy_param.set(adjncy_dummy);

        xadj_dummy.resize(2, 0);
        xadj_param.set(xadj_dummy);
    }

    // Resulting permutation and inverse permutation of the mesh
    idx_t perm[numCells];
    idx_t iperm[numCells];

    METIS_NodeND
    (
        &numCells,                       // num vertices in graph
        xadj_param.constCast().data(),   // indexing into adjncy
        adjncy_param.constCast().data(), // neighbour info
        NULL,
        NULL,
        perm,
        iperm
    );

    labelList newOrder(numCells);
    label cellInOrder = 0;
    for (int perm_idx = 0; perm_idx < numCells; perm_idx++)
    {
        newOrder[cellInOrder] = perm[perm_idx];
        cellInOrder++;
    }

    return newOrder;
}


// ************************************************************************* //
