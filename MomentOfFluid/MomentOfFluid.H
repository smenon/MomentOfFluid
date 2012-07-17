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
    MomentOfFluid

Description
    Moment-of-Fluid Interface Reconstruction as described in:

    H. T. Ahn, M. Shashkov,
    Multi-material interface reconstruction on generalized polyhedral meshes,
    Journal of Computational Physics,
    Volume 226, Issue 2, 1 October 2007, Pages 2096-2132,
    ISSN 0021-9991, 10.1016/j.jcp.2007.06.033.

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

SourceFiles
    MomentOfFluid.C

\*---------------------------------------------------------------------------*/

#ifndef MomentOfFluid_H
#define MomentOfFluid_H


#include "point.H"
#include "label.H"
#include "Tuple2.H"
#include "vector2D.H"
#include "FixedList.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;

/*---------------------------------------------------------------------------*\
                       Class MomentOfFluid Declaration
\*---------------------------------------------------------------------------*/

class MomentOfFluid
{
    // Private data

        //- Constant reference to mesh
        const polyMesh& mesh_;

        //- Hessian-normal plane definition
        typedef Tuple2<vector, scalar> hPlane;

        //- Triangle / Tetrahedron definition
        typedef FixedList<point, 3> Triangle;
        typedef FixedList<point, 4> Tetrahedron;

        //- Tet decomposition of original cell
        DynamicList<Tetrahedron> tetDecomp_;

        //- All intersection tets
        DynamicList<Tetrahedron> allTets_;

        //- Triangulated surfaces
        DynamicList<Triangle> allTris_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        MomentOfFluid(const MomentOfFluid&);

        //- Disallow default bitwise assignment
        void operator=(const MomentOfFluid&);

        //- Split and decompose
        void splitAndDecompose
        (
            const hPlane& clipPlane,
            const Tetrahedron& tet
        );

        // Extract triangles using plane info
        void extractTriangulation
        (
            const vector& xC,
            const hPlane& clipPlane,
            const Tetrahedron& tetra
        );

        //- Evaluate and return volume / centroid
        void getVolumeAndCentre(scalar& volume, vector& centre) const;

        // Class used during optimization
        class optInfo
        {
            // Reference to the parent class
            const MomentOfFluid& iRef_;

            // Cell index
            const label& cellIndex_;

            // Scalar fraction
            const scalar& fraction_;

            // Reference centroid
            const vector& refCentre_;

            // Resultant centroid
            vector& centre_;

            // Resultant distance
            scalar& distance_;

        public:

            // Constructor
            optInfo
            (
                const MomentOfFluid& iRef,
                const label& cellIndex,
                const scalar& fraction,
                const vector& refCentre,
                vector& centre,
                scalar& distance
            )
            :
                iRef_(iRef),
                cellIndex_(cellIndex),
                fraction_(fraction),
                refCentre_(refCentre),
                centre_(centre),
                distance_(distance)
            {}

            // Return reference to parent class
            const MomentOfFluid& reference() const
            {
                return iRef_;
            }

            const label& cellIndex() const
            {
                return cellIndex_;
            }

            const scalar& fraction() const
            {
                return fraction_;
            }

            const vector& refCentre() const
            {
                return refCentre_;
            }

            vector& centre()
            {
                return centre_;
            }

            scalar& distance()
            {
                return distance_;
            }
        };

        // Utility member to convert between csys
        inline vector sphericalToCartesian
        (
            const scalar& theta,
            const scalar& phi
        ) const
        {
            return
            (
                vector
                (
                    sin(theta) * cos(phi),
                    sin(theta) * sin(phi),
                    cos(theta)
                )
            );
        }

        // Integer sign
        inline label intSign(const scalar s) const
        {
            return (s >= 0.0)? 1: -1;
        }

        // Evaluate polynomial at location
        inline scalar evaluate
        (
            const FixedList<scalar, 4>& coeffs,
            const scalar x
        ) const
        {
            return
            (
                (coeffs[0] * x * x * x) +
                (coeffs[1] * x * x) +
                (coeffs[2] * x) +
                (coeffs[3])
            );
        }

        // Function evaluation routine
        scalar evaluateFunctional
        (
            const label& cellIndex,
            const scalar& fraction,
            const vector& refCentre,
            const vector2D& x,
            vector2D& fnGrad,
            vector& centre,
            scalar& distance
        );

        // Helper function for line-search
        scalar minimizeAlpha
        (
            const scalar endA,
            const scalar endB,
            const scalar alpha1,
            const scalar alpha2,
            const scalar f1,
            const scalar df1,
            const scalar f2,
            const scalar df2
        ) const;

        // Helper function for BFGS
        scalar lineSearch
        (
            const vector2D& x,
            const vector2D& dir,
            const scalar fInit,
            const scalar dfInit,
            const scalar alphaInit,
            const scalar rho,
            const scalar sigma,
            const scalar fMin,
            scalar& fAlpha,
            vector2D& gradAlpha,
            optInfo& data,
            label& flag,
            label& fnEvals
        );

        // Broyden-Fletcher-Goldfarb-Shanno (BFGS) algorithm
        scalar BFGS
        (
            vector2D& x,
            optInfo& data
        );

public:

    // Declare the name of the class and its debug switch

        TypeName("MomentOfFluid");

    // Constructors

        //- Construct from components
        MomentOfFluid(const polyMesh& mesh);


    // Destructor

        virtual ~MomentOfFluid();


    // Member Functions

        // Decompose original cell into tetrahedra
        void decomposeCell
        (
            const label& cellIndex
        );

        // Evaluate for intersections
        scalar evaluate
        (
            const Tuple2<vector, scalar>& plane,
            vector& centre
        );

        // Match specified volume fraction with supplied normal
        scalar matchFraction
        (
            const label& cellIndex,
            const scalar& fraction,
            const vector& normal,
            vector& centre,
            scalar& span,
            scalar* gdMin = NULL,
            scalar* gdMax = NULL
        );

        // Optimize for normal / centroid given a reference value
        void optimizeCentroid
        (
            const label& cellIndex,
            const scalar& fraction,
            const vector& refCentre,
            vector& normal,
            vector& centre
        );

        // Reconstruct the interface
        void constructInterface
        (
            const scalarField& fractions,
            const vectorField& refCentres
        );

        // Post-processing

            // Output trianglulated surface to VTK
            void outputSurface() const;

            // Output plane as VTK
            void outputPlane
            (
                const point& p,
                const vector& n,
                const label& cellIndex
            ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
