/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    tetIntersection

Description
    Tetrahedron-tetrahedron intersection

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

SourceFiles
    tetIntersectionI.H

\*---------------------------------------------------------------------------*/

#ifndef tetIntersection_H
#define tetIntersection_H

#include "MoF.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class tetIntersection Declaration
\*---------------------------------------------------------------------------*/

class tetIntersection
{
    // Private data

        //- Const reference to clipping tetrahedron
        const MoF::Tetrahedron& clipTet_;

        FixedList<MoF::hPlane, 4> clipPlanes_;

        //- Magnitude of clipping tetrahedron
        scalar clipTetMag_;

        //- Tetrahedra used as temporaries
        DynamicList<MoF::Tetrahedron> inside_;

        //- All intersection tets
        DynamicList<MoF::Tetrahedron> allTets_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        tetIntersection(const tetIntersection&);

        //- Disallow default bitwise assignment
        void operator=(const tetIntersection&);

        //- Compute clip-planes
        inline void computeClipPlanes();

public:

    // Constructors

        //- Construct from components
        inline tetIntersection(const MoF::Tetrahedron& clipTet);


    // Destructor

        inline ~tetIntersection();


    // Member Functions

        //- Return magnitude of clipping tetrahedron
        inline scalar clipTetMag() const;

        //- Evaluate for intersections against input tetrahedron
        inline bool evaluate(const MoF::Tetrahedron& subjectTet);

        //- Return intersections
        inline const DynamicList<MoF::Tetrahedron>& getIntersection() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tetIntersectionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
