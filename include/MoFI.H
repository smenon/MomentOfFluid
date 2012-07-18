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

Implemented by
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*---------------------------------------------------------------------------*/

namespace Foam
{

namespace MoF
{

// Decompose original cell into tetrahedra
//  - Optionally transform points to a local
//    coordinate system with the origin at xT.
void decomposeCell
(
    const polyMesh& mesh,
    const pointField& points,
    const label cellIndex,
    const point& xC,
    DynamicList<Tetrahedron>& tetDecomp,
    point xT
)
{
    Tetrahedron tmpTetra;

    // Clear list
    tetDecomp.clear();

    // Fetch references to connectivity
    const faceList& faces = mesh.faces();
    const cell& dCell = mesh.cells()[cellIndex];

    // Check for tetrahedral cell
    if (dCell.size() == 4)
    {
        // Insert points of cell
        const face& firstFace = faces[dCell[0]];
        const face& secondFace = faces[dCell[1]];

        // Fill first three points
        tmpTetra[0] = (points[firstFace[0]] - xT);
        tmpTetra[1] = (points[firstFace[1]] - xT);
        tmpTetra[2] = (points[firstFace[2]] - xT);

        // Fill isolated fourth point
        forAll(secondFace, pointI)
        {
            if
            (
                secondFace[pointI] != firstFace[0] &&
                secondFace[pointI] != firstFace[1] &&
                secondFace[pointI] != firstFace[2]
            )
            {
                tmpTetra[3] = (points[secondFace[pointI]] - xT);
                break;
            }
        }

        // Add tet to decomposition list
        tetDecomp.append(tmpTetra);
    }
    else
    {
        // Decompose using face-cell decomposition
        tmpTetra[3] = xC - xT;

        forAll(dCell, faceI)
        {
            const face& checkFace = faces[dCell[faceI]];

            // Optimize for triangle faces
            if (checkFace.size() == 3)
            {
                tmpTetra[0] = (points[checkFace[0]] - xT);
                tmpTetra[1] = (points[checkFace[1]] - xT);
                tmpTetra[2] = (points[checkFace[2]] - xT);

                // Add tet to decomposition list
                tetDecomp.append(tmpTetra);
            }
            else
            {
                // Pre-fill face centroid
                tmpTetra[2] = (checkFace.centre(points) - xT);

                forAll(checkFace, pI)
                {
                    tmpTetra[0] = (points[checkFace[pI]] - xT);
                    tmpTetra[1] = (points[checkFace.nextLabel(pI)] - xT);

                    // Add tet to decomposition list
                    tetDecomp.append(tmpTetra);
                }
            }
        }
    }
}


//- Evaluate and return volume / centroid
void getVolumeAndCentre
(
    const DynamicList<Tetrahedron>& tets,
    scalar& volume,
    vector& centre
)
{
    volume = 0.0;
    centre = vector::zero;

    forAll(tets, tetI)
    {
        const Tetrahedron& t = tets[tetI];

        // Calculate volume (no check for orientation)
        scalar tV =
        (
            Foam::mag
            (
                (1.0/6.0) *
                (
                    ((t[1] - t[0]) ^ (t[2] - t[0])) & (t[3] - t[0])
                )
            )
        );

        // Calculate centroid
        vector tC = (0.25 * (t[0] + t[1] + t[2] + t[3]));

        volume += tV;
        centre += (tV * tC);
    }

    centre /= volume + VSMALL;
}


// Split and decompose tetrahedron with supplied plane using
// the tetrahedron / half-space algorithm given in:
//
//   D.H. Eberly, 3D Game Engine Design: A Practical Approach to Real-time
//   Computer Graphics, Morgan Kaufmann, 2001.
//
//   Geometric Tools, LLC
//   Distributed under the Boost Software License, Version 1.0.
//   http://www.boost.org/LICENSE_1_0.txt
//
void splitAndDecompose
(
    const hPlane& clipPlane,
    const Tetrahedron& tet,
    DynamicList<Tetrahedron>& decompTets
)
{
    MoF::Tetrahedron tmpTetra;
    MoF::Tetrahedron tetra(tet);

    FixedList<scalar, 4> C;
    FixedList<label, 4> pos, neg, zero;
    label i = 0, nPos = 0, nNeg = 0, nZero = 0;

    for (i = 0; i < 4; ++i)
    {
        // Compute distance to plane
        C[i] = (tetra[i] & clipPlane.first()) - clipPlane.second();

        if (C[i] > 0.0)
        {
            pos[nPos++] = i;
        }
        else
        if (C[i] < 0.0)
        {
            neg[nNeg++] = i;
        }
        else
        {
            zero[nZero++] = i;
        }
    }

    if (nNeg == 0)
    {
        return;
    }

    if (nPos == 0)
    {
        decompTets.append(tetra);
        return;
    }

    // Tetrahedron is split by plane.  Determine how it is split and how to
    // decompose the negative-side portion into tetrahedra (6 cases).
    scalar w0, w1, invCDiff;
    vector intp[4];

    if (nPos == 3)
    {
        // +++-
        for (i = 0; i < nPos; ++i)
        {
            invCDiff = (1.0 / (C[pos[i]] - C[neg[0]]));

            w0 = -C[neg[0]] * invCDiff;
            w1 = +C[pos[i]] * invCDiff;

            tetra[pos[i]] = (w0 * tetra[pos[i]]) + (w1 * tetra[neg[0]]);
        }

        decompTets.append(tetra);
    }
    else
    if (nPos == 2)
    {
        if (nNeg == 2)
        {
            // ++--
            for (i = 0; i < nPos; ++i)
            {
                invCDiff = (1.0 / (C[pos[i]] - C[neg[0]]));

                w0 = -C[neg[0]] * invCDiff;
                w1 = +C[pos[i]] * invCDiff;

                intp[i] = (w0 * tetra[pos[i]]) + (w1 * tetra[neg[0]]);
            }

            for (i = 0; i < nNeg; ++i)
            {
                invCDiff = (1.0 / (C[pos[i]] - C[neg[1]]));

                w0 = -C[neg[1]] * invCDiff;
                w1 = +C[pos[i]] * invCDiff;

                intp[i+2] = (w0 * tetra[pos[i]]) + (w1 * tetra[neg[1]]);
            }

            tetra[pos[0]] = intp[2];
            tetra[pos[1]] = intp[1];

            decompTets.append(tetra);

            tmpTetra[0] = tetra[neg[1]];
            tmpTetra[1] = intp[3];
            tmpTetra[2] = intp[2];
            tmpTetra[3] = intp[1];

            decompTets.append(tmpTetra);

            tmpTetra[0] = tetra[neg[0]];
            tmpTetra[1] = intp[0];
            tmpTetra[2] = intp[1];
            tmpTetra[3] = intp[2];

            decompTets.append(tmpTetra);
        }
        else
        {
            // ++-0
            for (i = 0; i < nPos; ++i)
            {
                invCDiff = (1.0 / (C[pos[i]] - C[neg[0]]));

                w0 = -C[neg[0]] * invCDiff;
                w1 = +C[pos[i]] * invCDiff;

                tetra[pos[i]] = (w0 * tetra[pos[i]]) + (w1 * tetra[neg[0]]);
            }

            decompTets.append(tetra);
        }
    }
    else
    if (nPos == 1)
    {
        if (nNeg == 3)
        {
            // +---
            for (i = 0; i < nNeg; ++i)
            {
                invCDiff = (1.0 / (C[pos[0]] - C[neg[i]]));

                w0 = -C[neg[i]] * invCDiff;
                w1 = +C[pos[0]] * invCDiff;

                intp[i] = (w0 * tetra[pos[0]]) + (w1 * tetra[neg[i]]);
            }

            tetra[pos[0]] = intp[0];

            decompTets.append(tetra);

            tmpTetra[0] = intp[0];
            tmpTetra[1] = tetra[neg[1]];
            tmpTetra[2] = tetra[neg[2]];
            tmpTetra[3] = intp[1];

            decompTets.append(tmpTetra);

            tmpTetra[0] = tetra[neg[2]];
            tmpTetra[1] = intp[1];
            tmpTetra[2] = intp[2];
            tmpTetra[3] = intp[0];

            decompTets.append(tmpTetra);
        }
        else
        if (nNeg == 2)
        {
            // +--0
            for (i = 0; i < nNeg; ++i)
            {
                invCDiff = (1.0 / (C[pos[0]] - C[neg[i]]));

                w0 = -C[neg[i]] * invCDiff;
                w1 = +C[pos[0]] * invCDiff;

                intp[i] = (w0 * tetra[pos[0]]) + (w1 * tetra[neg[i]]);
            }

            tetra[pos[0]] = intp[0];

            decompTets.append(tetra);

            tmpTetra[0] = intp[1];
            tmpTetra[1] = tetra[zero[0]];
            tmpTetra[2] = tetra[neg[1]];
            tmpTetra[3] = intp[0];

            decompTets.append(tmpTetra);
        }
        else
        {
            // +-00
            invCDiff = (1.0 / (C[pos[0]] - C[neg[0]]));

            w0 = -C[neg[0]] * invCDiff;
            w1 = +C[pos[0]] * invCDiff;

            tetra[pos[0]] = (w0 * tetra[pos[0]]) + (w1 * tetra[neg[0]]);

            decompTets.append(tetra);
        }
    }
}


} // End namespace MoF


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
