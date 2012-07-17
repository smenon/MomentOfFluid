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
    testMomentOfFluid

Description
    Driver for testing MomentOfFluid intersection algorithms

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*---------------------------------------------------------------------------*/

#include "Time.H"
#include "fvMesh.H"
#include "argList.H"
#include "volFields.H"
#include "MomentOfFluid.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void unitTests(const polyMesh& mesh)
{
    vector refCentre;
    label testIndex = 12;
    vector testNormal(0.8427, 0.4815, 0.2408);

    // Normalize
    testNormal /= mag(testNormal);

    vector cellCentre = mesh.cellCentres()[testIndex];
    scalar cellVolume = mesh.cellVolumes()[testIndex];

    Info<< " Initial: " << nl
        << "  Normal: "
        << testNormal << nl
        << "  Cell Volume: " << cellVolume << nl
        << "  Cell Labels: "
        << mesh.cells()[testIndex].labels(mesh.faces()) << nl
        << "  Cell Points: "
        << mesh.cells()[testIndex].points(mesh.faces(), mesh.points()) << nl
        << endl;

    // Construct intersector
    MomentOfFluid mof(mesh);

    // Decompose original cell into tetrahedra
    mof.decomposeCell(testIndex);

    // Evaluate with specified plane
    scalar volume =
    (
        mof.evaluate
        (
            Tuple2<vector, scalar>
            (
                testNormal,
                0.0
            ),
            refCentre
        )
    );

    refCentre += cellCentre;

    scalar fraction = (volume / cellVolume);

//    // Match volume fraction with supplied normal
//    scalar span;
//
//    scalar d =
//    (
//        mof.matchFraction
//        (
//            testIndex,
//            testFraction,
//            testNormal,
//            centre,
//            span
//        )
//    );
//
//    scalar h = (0.01 * span);
//    scalar dMin = (d - h), dMax = (d + h);
//
//    mof.matchFraction
//    (
//        testIndex,
//        testFraction,
//        testNormal,
//        centre,
//        span,
//        &dMin,
//        &dMax
//    );

    vector normal = vector::zero, centre = vector::zero;

    mof.optimizeCentroid
    (
        testIndex,
        fraction,
        refCentre,
        normal,
        centre
    );

    Info<< " Final: " << nl
        << "  Normal: " << normal << nl
        << "  Centre: " << centre << nl
        << endl;

    mof.outputPlane(centre, normal, testIndex);
}

// Main program:

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

//    // Perform unit tests
//    unitTests(mesh);

    // Create fields
    volScalarField alpha
    (
        IOobject
        (
            "alpha1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    volVectorField refCentres
    (
        IOobject
        (
            "refCentres",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    // Construct intersector
    MomentOfFluid mof(mesh);

    // Compute surfaces and output
    mof.constructInterface
    (
        alpha.internalField(),
        refCentres.internalField()
    );

    // Output VTK file
    mof.outputSurface();
}
