/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam::COLAMDRenumber

Description
     Column approximate minimum degree renumbering

SourceFiles
    COLAMDRenumber.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_COLAMDRenumber_H
#define Foam_COLAMDRenumber_H

#include "renumberMethod.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class COLAMDRenumber Declaration
\*---------------------------------------------------------------------------*/

class COLAMDRenumber
:
    public renumberMethod
{
    // Private Member Functions

        //- No copy construct
        COLAMDRenumber(const COLAMDRenumber&) = delete;

        //- No copy assignment
        void operator=(const COLAMDRenumber&) = delete;


public:

    //- Runtime type information
    TypeName("COLAMD");


    // Constructors

        //- Construct given the renumber dictionary
        explicit COLAMDRenumber(const dictionary& dict);


    //- Destructor
    virtual ~COLAMDRenumber() = default;


    // Member Functions

        //- Return the order in which cells need to be visited
        //- (ie. from ordered back to original cell label).
        //  This is only defined for geometric renumberMethods.
        virtual labelList renumber(const pointField&) const
        {
            NotImplemented;
            return labelList();
        }

        //- Return the order in which cells need to be visited
        //- (ie. from ordered back to original cell label).
        //  Use the mesh connectivity (if needed)
        virtual labelList renumber
        (
            const polyMesh& mesh,
            const pointField& cc
        ) const;

        //- Return the order in which cells need to be visited
        //- (ie. from ordered back to original cell label).
        //  Connectivity in losort addressing (= neighbour + offsets into
        //  neighbour)
        virtual labelList renumber
        (
            const labelList& cellCells,
            const labelList& offsets,
            const pointField& cc
        ) const;

        //- Return the order in which cells need to be visited
        //- (ie. from ordered back to original cell label).
        virtual labelList renumber
        (
            const CompactListList<label>& cellCells,
            const pointField& cellCentres
        ) const;

        //- Return the order in which cells need to be visited
        //- (ie. from ordered back to original cell label).
        //  The connectivity is equal to mesh.cellCells() except
        //  - the connections are across coupled patches
        virtual labelList renumber
        (
            const labelListList& cellCells,
            const pointField& cellCentres
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
