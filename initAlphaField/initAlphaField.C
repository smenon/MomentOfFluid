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
    // Set reference to target points / cells
    const cellList& tgtCells = meshTarget.cells();
    const pointField& tgtPoints = meshTarget.points();

    // Fetch references to centres
    const pointField& srcCentres = meshSource.cellCentres();
    const pointField& tgtCentres = meshTarget.cellCentres();

    // Set reference to source cellCells
    const labelListList& cc = meshSource.cellCells();

    // Construct and fetch source cell-search tree
    const indexedOctree<treeDataCell>& tree = meshSource.cellTree();

    // Tet decomposition of cells
    DynamicList<MoF::Tetrahedron> srcDecomp(10);
    DynamicList<MoF::Tetrahedron> tgtDecomp(10);

    forAll(tgtCells, cellI)
    {
        const cell& tgtCell = tgtCells[cellI];

        // Find nearest source cell
        label nearest = tree.findInside(tgtCentres[cellI]);

        if (nearest == -1)
        {
            const labelList pLabels(tgtCell.labels(meshTarget.faces()));

            forAll(pLabels, pI)
            {
                if ((nearest = tree.findInside(tgtPoints[pLabels[pI]])) > -1)
                {
                    break;
                }
            }
        }

        // Assume outside and skip
        if (nearest == -1)
        {
            alpha.internalField()[cellI] = 0.0;
            refCentres.internalField()[cellI] = vector::zero;

            continue;
        }

        // Found a candidate, commence search
        scalar tgtAlpha = 0.0;
        vector tgtRefCentre = vector::zero;

        // Fetch volume
        scalar tgtVolume = meshTarget.cellVolumes()[cellI];

        // Decompose target cell, if necessary
        MoF::decomposeCell
        (
            meshTarget,
            meshTarget.points(),
            cellI,
            tgtCentres[cellI],
            tgtDecomp
        );

        // Initialize target intersectors
        PtrList<tetIntersection> tgtInt(tgtDecomp.size());

        forAll(tgtInt, intI)
        {
            tgtInt.set(intI, new tetIntersection(tgtDecomp[intI]));
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

                    // Decompose source cell, if necessary.
                    MoF::decomposeCell
                    (
                        meshSource,
                        meshSource.points(),
                        checkEntity,
                        srcCentres[checkEntity],
                        srcDecomp
                    );

                    bool anyIntersects = false;

                    forAll(srcDecomp, tetI)
                    {
                        const MoF::Tetrahedron& tetraI = srcDecomp[tetI];

                        forAll(tgtDecomp, tetJ)
                        {
                            // Intersect source / target tets
                            tetIntersection& tJ = tgtInt[tetJ];

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
                                tgtAlpha += volume;
                                tgtRefCentre += (volume * centre);

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

        // Normalize
        tgtRefCentre /= tgtAlpha + VSMALL;
        tgtAlpha /= tgtVolume + VSMALL;

        // Fix over-fill
        if (tgtAlpha > 1.0)
        {
            tgtAlpha = 1.0;
        }

        // Set accumulate result
        alpha.internalField()[cellI] = tgtAlpha;
        refCentres.internalField()[cellI] = tgtRefCentre;
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
