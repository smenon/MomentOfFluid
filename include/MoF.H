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

Namespace
    MoF

Description
    Utility functions for Moment of Fluid algorithms

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef MoF_H
#define MoF_H

#include "point.H"
#include "label.H"
#include "Tuple2.H"
#include "polyMesh.H"
#include "FixedList.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Namespace MoF Declaration
\*---------------------------------------------------------------------------*/

namespace MoF
{
    //- Hessian-normal plane definition
    typedef Tuple2<vector, scalar> hPlane;

    //- Triangle / Tetrahedron definition
    typedef FixedList<point, 3> Triangle;
    typedef FixedList<point, 4> Tetrahedron;

    //- Decompose cell into tetrahedra
    void decomposeCell
    (
        const polyMesh& mesh,
        const pointField& points,
        const label cellIndex,
        const point& xC,
        DynamicList<Tetrahedron>& tetDecomp,
        point xT = vector::zero
    );

    //- Evaluate and return volume / centroid
    void getVolumeAndCentre
    (
        const DynamicList<Tetrahedron>& tets,
        scalar& volume,
        vector& centre
    );

    //- Split and decompose
    void splitAndDecompose
    (
        const hPlane& clipPlane,
        const Tetrahedron& tet,
        DynamicList<Tetrahedron>& decompTets
    );

} // End namespace MoF


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "MoFI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
