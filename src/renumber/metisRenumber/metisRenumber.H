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
    Foam::metisRenumber

Description
    Use Multilevel Nested Dissection for fill reducing sparse matrix ordering

SourceFiles
    metisRenumber.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_metisRenumber_H
#define Foam_metisRenumber_H

#include "renumberMethod.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class metisRenumber Declaration
\*---------------------------------------------------------------------------*/

class metisRenumber
:
    public renumberMethod
{
    // Private Data

        const dictionary& coeffsDict_;


    // Private Member Functions

        //- No copy construct
        metisRenumber(const metisRenumber&) = delete;

        //- No copy assignment
        void operator=(const metisRenumber&) = delete;


public:

    //- Runtime type information
    TypeName("metis");


    // Constructors

        //- Construct given the renumber dictionary
        explicit metisRenumber(const dictionary& dict);


    //- Destructor
    virtual ~metisRenumber() = default;


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
        virtual labelList renumber
        (
            const CompactListList<label>& cellCells,
            const pointField& cellCentres
        ) const
        {
            NotImplemented;
            return labelList();
        }

        //- Return the order in which cells need to be visited
        //- (ie. from ordered back to original cell label).
        virtual labelList renumber
        (
            const labelListList& cellCells,
            const pointField& cellCentres
        ) const
        {
            NotImplemented;
            return labelList();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
