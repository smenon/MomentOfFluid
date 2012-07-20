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

Application
    initAlphaField

Description
    Initialize fields for Moment-Of-Fluid interfaces

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*---------------------------------------------------------------------------*/

#include "Time.H"
#include "fvCFD.H"
#include "argList.H"
#include "treeDataCell.H"
#include "pointIndexHit.H"
#include "indexedOctree.H"
#include "tetIntersection.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Calculate and populate fields
void initAlphaField
(
    const fvMesh& meshSource,
    const fvMesh& meshTarget,
    volScalarField& alpha,
    volVectorField& refCentres
)
{
    // Initialize fields
    scalarField& aiF = alpha.internalField();
    vectorField& rCiF = refCentres.internalField();

    aiF = 0.0;
    rCiF = vector::zero;

    // Set reference to target points / cells
    const cellList& tgtCells = meshTarget.cells();
    const pointField& tgtPoints = meshTarget.points();

    // Set reference to source points / cells
    const cellList& srcCells = meshSource.cells();
    const pointField& srcPoints = meshSource.points();

    // Fetch references to centres / volumes
    const pointField& srcCentres = meshSource.cellCentres();
    const pointField& tgtCentres = meshTarget.cellCentres();

    const scalarField& srcVolumes = meshSource.cellVolumes();
    const scalarField& tgtVolumes = meshTarget.cellVolumes();

    // Set reference to target cellCells
    const labelListList& cc = meshTarget.cellCells();

    // Construct and fetch target cell-search tree
    const indexedOctree<treeDataCell>& tree = meshTarget.cellTree();

    // Tet decomposition of cells
    DynamicList<MoF::Tetrahedron> srcDecomp(10);
    DynamicList<MoF::Tetrahedron> tgtDecomp(10);

    forAll(srcCells, cellI)
    {
        // Find nearest target cell
        label nearest = tree.findInside(srcCentres[cellI]);

        if (nearest == -1)
        {
            FatalErrorIn
            (
                "\n\n"
                "void initAlphaField\n"
                "(\n"
                "    const fvMesh& meshSource,\n"
                "    const fvMesh& meshTarget,\n"
                "    volScalarField& alpha,\n"
                "    volVectorField& refCentres\n"
                ")\n"
            )
                << " Suitable target cell was not found." << nl
                << " Index: " << cellI << nl
                << abort(FatalError);
        }

        // Fetch volume
        scalar volAlpha = 0.0;
        scalar srcVolume = srcVolumes[cellI];

        // Decompose source cell, if necessary
        MoF::decomposeCell
        (
            meshSource,
            srcPoints,
            cellI,
            srcCentres[cellI],
            srcDecomp
        );

        // Initialize source intersectors
        PtrList<tetIntersection> srcInt(srcDecomp.size());

        forAll(srcInt, intI)
        {
            srcInt.set(intI, new tetIntersection(srcDecomp[intI]));
        }

        bool changed;
        label nAttempts = 0;
        labelHashSet checked, skipped;

        // Loop and add intersections until nothing changes
        do
        {
            // Reset flag
            changed = false;

            // Fetch the set of candidates
            labelList checkList;

            if (nAttempts == 0)
            {
                checkList = labelList(1, nearest);
            }
            else
            {
                checkList = checked.toc();
            }

            forAll(checkList, indexI)
            {
                labelList checkEntities;

                if (nAttempts == 0)
                {
                    checkEntities = labelList(1, checkList[indexI]);
                }
                else
                {
                    checkEntities = cc[checkList[indexI]];
                }

                forAll(checkEntities, entityI)
                {
                    label checkEntity = checkEntities[entityI];

                    // Skip if this is already
                    // on the checked / skipped list
                    if
                    (
                        (checked.found(checkEntity)) ||
                        (skipped.found(checkEntity))
                    )
                    {
                        continue;
                    }

                    // Decompose target cell, if necessary.
                    MoF::decomposeCell
                    (
                        meshTarget,
                        tgtPoints,
                        checkEntity,
                        tgtCentres[checkEntity],
                        tgtDecomp
                    );

                    bool anyIntersects = false;

                    forAll(tgtDecomp, tetI)
                    {
                        const MoF::Tetrahedron& tetraI = tgtDecomp[tetI];

                        forAll(srcDecomp, tetJ)
                        {
                            // Intersect source / target tets
                            tetIntersection& tJ = srcInt[tetJ];

                            bool intersect = tJ.evaluate(tetraI);

                            if (intersect)
                            {
                                scalar volume = 0.0;
                                vector centre = vector::zero;

                                // Get volume / centroid
                                MoF::getVolumeAndCentre
                                (
                                    tJ.getIntersection(),
                                    volume,
                                    centre
                                );

                                // Accumulate result
                                volAlpha += volume;
                                aiF[checkEntity] += volume;
                                rCiF[checkEntity] += (volume * centre);

                                anyIntersects = true;
                            }
                        }
                    }

                    if (anyIntersects)
                    {
                        if (!checked.found(checkEntity))
                        {
                            checked.insert(checkEntity);
                        }

                        changed = true;
                    }
                    else
                    if (nAttempts == 0)
                    {
                        FatalErrorIn
                        (
                            "\n\n"
                            "void initAlphaField\n"
                            "(\n"
                            "    const fvMesh& meshSource,\n"
                            "    const fvMesh& meshTarget,\n"
                            "    volScalarField& alpha,\n"
                            "    volVectorField& refCentres\n"
                            ")\n"
                        )
                            << " First intersection was not found." << nl
                            << " Index: " << cellI << nl
                            << " nearest: " << nearest << nl
                            << " checkEntity: " << checkEntity << nl
                            << " srcDecomp: " << srcDecomp << nl
                            << " tgtDecomp: " << tgtDecomp << nl
                            << abort(FatalError);
                    }
                    else
                    {
                        // Add to the skipped list
                        if (!skipped.found(checkEntity))
                        {
                            skipped.insert(checkEntity);
                        }
                    }
                }
            }

            nAttempts++;

            // Break out if we're taking too long
            if (nAttempts > 20)
            {
                break;
            }

        } while (changed);

        // Check if volume was completely enclosed
        scalar error = Foam::mag(1.0 - (volAlpha / srcVolume));

        if (error > 1e-10)
        {
            Info<< " Cell: " << cellI << nl
                << "  Error: " << error << nl
                << "  volAlpha: " << volAlpha << nl
                << "  srcVolume: " << srcVolume << nl
                << abort(FatalError);
        }
    }

    // Normalize
    forAll(tgtCells, cellI)
    {
        scalar& alpha = aiF[cellI];
        vector& refCentre = rCiF[cellI];
        scalar tgtVolume = tgtVolumes[cellI];

        if (Foam::mag(alpha) > 0.0)
        {
            refCentre /= alpha;
            alpha /= tgtVolume;
        }
    }
}


// Main program:
int main(int argc, char *argv[])
{
#   include "setRoots.H"
#   include "createTimes.H"
#   include "setTimeIndex.H"

    runTimeSource.setTime(sourceTimes[sourceTimeIndex], sourceTimeIndex);
    runTimeTarget.setTime(sourceTimes[sourceTimeIndex], sourceTimeIndex);

    Info<< "\nSource time: " << runTimeSource.value()
        << "\nTarget time: " << runTimeTarget.value()
        << endl;

    Info<< "Create meshes\n" << endl;

    fvMesh meshSource
    (
        IOobject
        (
            fvMesh::defaultRegion,
            runTimeSource.timeName(),
            runTimeSource
        )
    );

    fvMesh meshTarget
    (
        IOobject
        (
            fvMesh::defaultRegion,
            runTimeTarget.timeName(),
            runTimeTarget
        )
    );

    // Create fields
    volScalarField alpha
    (
        IOobject
        (
            "alpha1",
            runTimeTarget.timeName(),
            meshTarget,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        meshTarget,
        dimensionedScalar("alpha", dimless, 0.0)
    );

    volVectorField refCentres
    (
        IOobject
        (
            "refCentres",
            runTimeTarget.timeName(),
            meshTarget,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        meshTarget,
        dimensionedVector("refCentre", dimless, vector::zero)
    );

    // Calculate and populate fields
    initAlphaField
    (
        meshSource,
        meshTarget,
        alpha,
        refCentres
    );

    // Write fields
    alpha.write();
    refCentres.write();

    return 0;
}


// ************************************************************************* //
