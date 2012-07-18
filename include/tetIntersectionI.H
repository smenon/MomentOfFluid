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

\*---------------------------------------------------------------------------*/

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Compute clip-planes
inline void tetIntersection::computeClipPlanes()
{
    // Define edge vectors
    vector edge10 = clipTet_[1] - clipTet_[0];
    vector edge20 = clipTet_[2] - clipTet_[0];
    vector edge30 = clipTet_[3] - clipTet_[0];
    vector edge21 = clipTet_[2] - clipTet_[1];
    vector edge31 = clipTet_[3] - clipTet_[1];

    // Cross-products
    clipPlanes_[0].first() = (edge20 ^ edge10);
    clipPlanes_[1].first() = (edge10 ^ edge30);
    clipPlanes_[2].first() = (edge30 ^ edge20);
    clipPlanes_[3].first() = (edge21 ^ edge31);

    // Normalize
    clipPlanes_[0].first() /= mag(clipPlanes_[0].first()) + VSMALL;
    clipPlanes_[1].first() /= mag(clipPlanes_[1].first()) + VSMALL;
    clipPlanes_[2].first() /= mag(clipPlanes_[2].first()) + VSMALL;
    clipPlanes_[3].first() /= mag(clipPlanes_[3].first()) + VSMALL;

    // Compute magnitude of clipping tetrahedron
    clipTetMag_ = (1.0 / 6.0) * (edge10 & clipPlanes_[3].first());

    if (clipTetMag_ < 0.0)
    {
        // Reverse normal directions
        clipPlanes_[0].first() = -clipPlanes_[0].first();
        clipPlanes_[1].first() = -clipPlanes_[1].first();
        clipPlanes_[2].first() = -clipPlanes_[2].first();
        clipPlanes_[3].first() = -clipPlanes_[3].first();

        // Reverse sign
        clipTetMag_ = mag(clipTetMag_);
    }

    // Determine plane constants
    clipPlanes_[0].second() = (clipTet_[0] & clipPlanes_[0].first());
    clipPlanes_[1].second() = (clipTet_[1] & clipPlanes_[1].first());
    clipPlanes_[2].second() = (clipTet_[2] & clipPlanes_[2].first());
    clipPlanes_[3].second() = (clipTet_[3] & clipPlanes_[3].first());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline tetIntersection::tetIntersection(const FixedList<point, 4>& clipTet)
:
    clipTet_(clipTet),
    clipTetMag_(0.0),
    inside_(10),
    allTets_(10)
{
    // Pre-compute clipping planes
    computeClipPlanes();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

inline tetIntersection::~tetIntersection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return magnitude of clipping tetrahedron
inline scalar tetIntersection::clipTetMag() const
{
    return clipTetMag_;
}


// Evaluate for intersections
inline bool tetIntersection::evaluate(const FixedList<point, 4>& subjectTet)
{
    // Clear lists
    inside_.clear();
    allTets_.clear();

    // Add initial tetrahedron to list
    allTets_.append(subjectTet);

    // Clip second tetrahedron against planes of clipping tetrahedron
    for (label i = 0; i < 4; i++)
    {
        forAll(allTets_, tetI)
        {
            MoF::splitAndDecompose
            (
                clipPlanes_[i],
                allTets_[tetI],
                inside_
            );
        }

        // Prep for next clipping plane
        allTets_ = inside_;
        inside_.clear();
    }

    return (allTets_.size() > 0);
}


// Return intersections
inline const DynamicList<FixedList<point, 4> >&
tetIntersection::getIntersection() const
{
    return allTets_;
}


}

// ************************************************************************* //
