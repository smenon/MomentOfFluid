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
    Moment-of-Fluid Interface Reconstruction algorithm for general polyhedra

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*---------------------------------------------------------------------------*/

#include "Time.H"
#include "OFstream.H"
#include "tensor2D.H"

#include "MomentOfFluid.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(MomentOfFluid, 0);

// Defined at file scope
static scalar eps_ = 2.2204e-16;
static scalar sqrteps_ = 1.4901e-08;
static scalar cbrteps_ = 6.0554e-06;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Extract triangles using plane info
//  - Modified version of splitAndDecompose
void MomentOfFluid::extractTriangulation
(
    const vector& xC,
    const MoF::hPlane& clipPlane,
    const MoF::Tetrahedron& tetra
)
{
    MoF::Triangle tmpTri;

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

    if (nNeg == 0 || nPos == 0)
    {
        return;
    }

    // Tetrahedron is split by plane.
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

            tmpTri[i] = (w0 * tetra[pos[i]]) + (w1 * tetra[neg[0]]) + xC;
        }

        allTris_.append(tmpTri);
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

            tmpTri[0] = intp[3] + xC;
            tmpTri[1] = intp[2] + xC;
            tmpTri[2] = intp[1] + xC;

            allTris_.append(tmpTri);

            tmpTri[0] = intp[0] + xC;
            tmpTri[1] = intp[1] + xC;
            tmpTri[2] = intp[2] + xC;

            allTris_.append(tmpTri);
        }
        else
        {
            // ++-0
            for (i = 0; i < nPos; ++i)
            {
                invCDiff = (1.0 / (C[pos[i]] - C[neg[0]]));

                w0 = -C[neg[0]] * invCDiff;
                w1 = +C[pos[i]] * invCDiff;

                tmpTri[i] = (w0 * tetra[pos[i]]) + (w1 * tetra[neg[0]]) + xC;
            }

            tmpTri[2] = tetra[zero[0]] + xC;

            allTris_.append(tmpTri);
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

                tmpTri[i] = (w0 * tetra[pos[0]]) + (w1 * tetra[neg[i]]) + xC;
            }

            allTris_.append(tmpTri);
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

                tmpTri[i] = (w0 * tetra[pos[0]]) + (w1 * tetra[neg[i]]) + xC;
            }

            tmpTri[2] = tetra[zero[0]] + xC;

            allTris_.append(tmpTri);
        }
        else
        {
            // +-00
            invCDiff = (1.0 / (C[pos[0]] - C[neg[0]]));

            w0 = -C[neg[0]] * invCDiff;
            w1 = +C[pos[0]] * invCDiff;

            tmpTri[0] = (w0 * tetra[pos[0]]) + (w1 * tetra[neg[0]]) + xC;
            tmpTri[1] = tetra[zero[0]] + xC;
            tmpTri[2] = tetra[zero[1]] + xC;

            allTris_.append(tmpTri);
        }
    }
}


// Evaluate for intersections
scalar MomentOfFluid::evaluate
(
    const Tuple2<vector, scalar>& plane,
    vector& centre
)
{
    scalar volume = 0.0;

    // Clear list
    allTets_.clear();

    // Clip tetrahedra against the plane
    forAll(tetDecomp_, tetI)
    {
        MoF::splitAndDecompose
        (
            plane,
            tetDecomp_[tetI],
            allTets_
        );
    }

    // Compute quantities
    MoF::getVolumeAndCentre(allTets_, volume, centre);

    return volume;
}


// Function evaluation routine
scalar MomentOfFluid::evaluateFunctional
(
    const label& cellIndex,
    const scalar& fraction,
    const vector& refCentre,
    const vector2D& x,
    vector2D& fnGrad,
    vector& centre,
    scalar& distance
)
{
    scalar fnVal = 0.0, span = 0.0;

    // Recover the normal from inputs
    vector normal = sphericalToCartesian(x[0], x[1]);

    // Match fraction with supplied normal
    distance =
    (
        matchFraction
        (
            cellIndex,
            fraction,
            normal,
            centre,
            span
        )
    );

    // Evaluate functional
    fnVal = Foam::mag(refCentre - centre);

    // Optimize search for gradient steps
    //scalar hg = (0.01 * span);
    //scalar dMin = (distance - hg), dMax = (distance + hg);

    // Evaluate a finite-difference gradient
    bool central = true;

    // Specify step-size for gradient
    scalar h = central ? cbrteps_ : sqrteps_;

    scalar fVal[2];
    vector fC[2], fNorm[2];

    fNorm[0] = sphericalToCartesian(x[0] + h, x[1]);
    fNorm[1] = sphericalToCartesian(x[0], x[1] + h);

    // Evaluate forward difference
    for (label i = 0; i < 2; i++)
    {
        matchFraction
        (
            cellIndex,
            fraction,
            fNorm[i],
            fC[i],
            span
        );

        fVal[i] = Foam::mag(refCentre - fC[i]);
    }

    // Optionally evaluate central differences
    if (central)
    {
        scalar rVal[2];
        vector rC[2], rNorm[2];

        rNorm[0] = sphericalToCartesian(x[0] - h, x[1]);
        rNorm[1] = sphericalToCartesian(x[0], x[1] - h);

        for (label i = 0; i < 2; i++)
        {
            matchFraction
            (
                cellIndex,
                fraction,
                rNorm[i],
                rC[i],
                span
            );

            rVal[i] = Foam::mag(refCentre - rC[i]);
        }

        fnGrad[0] = (fVal[0] - rVal[0]) / (2 * h);
        fnGrad[1] = (fVal[1] - rVal[1]) / (2 * h);
    }
    else
    {
        fnGrad[0] = (fVal[0] - fnVal) / h;
        fnGrad[1] = (fVal[1] - fnVal) / h;
    }

    return fnVal;
}


// Match specified volume fraction with supplied normal
//  - Supplied normal is expected to be normalized
//  - Return distance corresponding to matched fraction
//  - Optionally use supplied guesses to improve convergence
scalar MomentOfFluid::matchFraction
(
    const label& cellIndex,
    const scalar& fraction,
    const vector& normal,
    vector& centre,
    scalar& span,
    scalar* gdMin,
    scalar* gdMax
)
{
    // Fetch cell volume / centroid
    const vector& xC = mesh_.cellCentres()[cellIndex];
    const scalar& volume = mesh_.cellVolumes()[cellIndex];

    // Compute limits of signed distance along normal
    label fEvals = 0;
    scalar fdMin = 0.0, fdMax = 1.0;
    scalar dMin = GREAT, dMax = -GREAT;

    if (gdMin == NULL && gdMax == NULL)
    {
        forAll(tetDecomp_, tetI)
        {
            const MoF::Tetrahedron& tet = tetDecomp_[tetI];

            forAll(tet, pointI)
            {
                scalar project = (tet[pointI] & normal);

                dMin = Foam::min(dMin, project);
                dMax = Foam::max(dMax, project);
            }
        }
    }
    else
    {
        // Use supplied guesses
        dMin = *gdMin;
        dMax = *gdMax;

        // Evaluate function at guesses
        fdMin = (evaluate(MoF::hPlane(normal, dMin), centre) / volume);
        fdMax = (evaluate(MoF::hPlane(normal, dMax), centre) / volume);

        fEvals += 2;
    }

    // Specify span
    span = (dMax - dMin);

    // Secant-method for volume-matching
    scalar error, tol = 1e-10;
    label iter = 0, maxIter = 50;

    scalar fd, d, deltaf;
    scalar d1 = dMax, d2 = dMin;
    scalar f1 = fdMax - fraction, f2 = fdMin - fraction;

    bool useBisection = false;

    do
    {
        deltaf = (f1 - f2);

        if (!useBisection && mag(deltaf) < VSMALL)
        {
            // Fallback to bisection
            useBisection = true;

            iter = 0;
            d1 = dMax; d2 = dMin;
            f1 = fdMax - fraction; f2 = fdMin - fraction;
        }

        if (useBisection)
        {
            // Update by bisection method
            d = 0.5 * (d1 + d2);
        }
        else
        {
            // Update distance by secant method
            d = d1 - f1 * (d1 - d2) / deltaf;
        }

        // Compute functional
        fd = (evaluate(MoF::hPlane(normal, d), centre) / volume) - fraction;

        // Compute error
        error = Foam::mag(fd);

        // Update for next recurrence
        if (useBisection)
        {
            if (intSign(fd) == intSign(f1))
            {
                d1 = d;
                f1 = fd;
            }
            else
            if (intSign(fd) == intSign(f2))
            {
                d2 = d;
                f2 = fd;
            }
            else
            {
                FatalErrorIn("void MomentOfFluid::matchFraction()")
                    << " Bisection not feasible." << nl
                    << "   cellIndex: " << cellIndex << nl
                    << "   fraction: " << fraction << nl
                    << "   fd: " << fd << nl
                    << "   f1: " << f1 << nl
                    << "   f2: " << f2 << nl
                    << abort(FatalError);
            }
        }
        else
        {
            // Update for secant method
            d2 = d1; f2 = f1;
            d1 = d; f1 = fd;
        }

        iter++;
        fEvals++;

    } while (iter < maxIter && error > tol);

    if (iter == maxIter)
    {
        InfoIn("void MomentOfFluid::matchFraction()")
            << nl << " Max iterations reached. "
            << nl << "   cellIndex: " << cellIndex
            << nl << "   fraction: " << fraction
            << nl << "   error: " << error
            << endl;
    }

    // Add centroid to result
    centre += xC;

    if (debug > 1)
    {
        Info<< " Iter: " << iter
            << " fEvals: " << fEvals << nl
            << "   Distance: " << d << nl
            << "   refCentre: " << centre << nl
            << "   Centre: " << xC << nl
            << "   Normal:" << normal << nl
            << "   Error: " << error
            << "   Span: " << span << nl
            << endl;
    }

    return d;
}


// Optimize for normal / centroid given a reference value
void MomentOfFluid::optimizeCentroid
(
    const label& cellIndex,
    const scalar& fraction,
    const vector& refCentre,
    vector& normal,
    vector& centre
)
{
    const vector& xC = mesh_.cellCentres()[cellIndex];

    // Decompose cell, transforming to cell centroid
    MoF::decomposeCell
    (
        mesh_,
        mesh_.points(),
        cellIndex,
        xC,
        tetDecomp_,
        xC
    );

    scalar distance = 0.0;

    // Make an initial guess for the normal
    vector iNormal = (xC - refCentre);

    iNormal /= mag(iNormal) + VSMALL;

    // Convert components to spherical coordinates
    scalar theta = acos(iNormal.z());
    scalar phi = atan2(iNormal.y(), iNormal.x());

    // Prepare inputs to BFGS
    vector2D x(theta, phi);

    // Prepare data
    optInfo data
    (
        *this,
        cellIndex,
        fraction,
        refCentre,
        centre,
        distance
    );

    if (debug)
    {
        Info<< " Initial: " << iNormal << nl
            << "   Theta: " << theta << nl
            << "   Phi: " << phi << endl;
    }

    // Call the bfgs algorithm
    scalar fnVal = BFGS(x, data);

    // Update result
    normal = sphericalToCartesian(x[0], x[1]);

    if (debug)
    {
        Info<< " Final: " << nl
            << "  Functional: " << fnVal << nl
            << "  Normal: " << normal << nl
            << "  Centre: " << centre << nl
            << "  Distance: " << distance << nl
            << endl;
    }

    // Output triangulation
    if (debug)
    {
        forAll(tetDecomp_, tetI)
        {
            extractTriangulation
            (
                xC,
                MoF::hPlane(normal, distance),
                tetDecomp_[tetI]
            );
        }
    }
}


// Helper function for line-search
scalar MomentOfFluid::minimizeAlpha
(
    const scalar endA,
    const scalar endB,
    const scalar alpha1,
    const scalar alpha2,
    const scalar f1,
    const scalar df1,
    const scalar f2,
    const scalar df2
) const
{
    scalar alpha = 0.0, dAlpha = (alpha2 - alpha1);

    // Specify cubic polynomial coefficients for interpolation
    FixedList<scalar, 2> roots;
    FixedList<scalar, 4> coeffs;

    coeffs[3] = f1;
    coeffs[2] = dAlpha * df1;
    coeffs[1] = 3.0 * (f2 - f1) - (2.0 * df1 + df2) * dAlpha;
    coeffs[0] = (df1 + df2) * dAlpha - 2.0 * (f2 - f1);

    // Specify bounds
    scalar lb = (endA - alpha1) / (alpha2 - alpha1);
    scalar ub = (endB - alpha1) / (alpha2 - alpha1);

    // Swap bounds if necessary
    if (lb > ub)
    {
        Foam::Swap(lb, ub);
    }

    // Evaluate at bounds
    scalar lbp = evaluate(coeffs, lb);
    scalar ubp = evaluate(coeffs, ub);

    scalar z = (lbp < ubp) ? lb : ub;
    scalar pmin = Foam::min(lbp, ubp);

    // Find roots with specified coefficients
    scalar a = 3.0 * coeffs[0];
    scalar b = 2.0 * coeffs[1];
    scalar c = 1.0 * coeffs[2];

    // Determine discriminant
    scalar delta = (b * b) - (4.0 * a * c);

    // Check roots if they're real
    if (delta >= 0.0)
    {
        if (mag(a) < VSMALL)
        {
            roots[0] = roots[1] = (-c / b);
        }
        else
        {
            roots[0] = (-b + Foam::sqrt(delta)) / (2.0 * a);
            roots[1] = (-b - Foam::sqrt(delta)) / (2.0 * a);
        }

        if (lb <= roots[1] && roots[1] <= ub)
        {
            scalar r1p = evaluate(coeffs, roots[1]);

            z = (r1p < pmin) ? roots[1] : z;
            pmin = Foam::min(pmin, r1p);
        }

        if (lb <= roots[0] && roots[0] <= ub)
        {
            scalar r0p = evaluate(coeffs, roots[0]);

            z = (r0p < pmin) ? roots[0] : z;
            pmin = Foam::min(pmin, r0p);
        }
    }

    // Evaluate alpha and return
    alpha = alpha1 + z * (alpha2 - alpha1);

    return alpha;
}


// Helper function for BFGS
scalar MomentOfFluid::lineSearch
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
)
{
    scalar alpha = alphaInit;

    // Step 1: Look for an acceptable bracket
    fAlpha = fInit;

    scalar fPrev, dfPrev;
    scalar dfAlpha = dfInit;
    scalar alphaPrev = 0.0, alphaMax = (fMin - fInit) / (rho * dfInit);

    // Brackets, functions and derivatives
    scalar a = 0.0, b = 0.0;
    scalar fa = 0.0, fb = 0.0;
    scalar dfa = 0.0, dfb = 0.0;

    // Parameters
    scalar t1 = 9.0, t2 = Foam::min(0.1, sigma), t3 = 0.5;

    label iter = 0, maxIter = 100;

    while (iter < maxIter)
    {
        fPrev = fAlpha;
        dfPrev = dfAlpha;

        // Evaluate function and gradient
        fAlpha =
        (
            evaluateFunctional
            (
                data.cellIndex(),
                data.fraction(),
                data.refCentre(),
                (x + (alpha * dir)),
                gradAlpha,
                data.centre(),
                data.distance()
            )
        );

        dfAlpha = (gradAlpha & dir);

        iter++;
        fnEvals++;

        if (fAlpha < fMin)
        {
            flag = 1;
            break;
        }

        if (fAlpha > (fInit + alpha * rho * dfInit) || fAlpha >= fPrev)
        {
            a = alphaPrev; b = alpha;
            fa = fPrev; dfa = dfPrev;
            fb = fAlpha; dfb = dfAlpha;
            flag = 2;
            break;
        }

        if (mag(dfAlpha) <= -sigma * dfInit)
        {
            flag = 0;
            break;
        }

        if (dfAlpha >= 0.0)
        {
            a = alpha; b = alphaPrev;
            fa = fAlpha; dfa = dfAlpha;
            fb = fPrev; dfb = dfPrev;
            flag = 2;
            break;
        }

        // Update alpha
        if ((2.0 * alpha) - alphaPrev < alphaMax)
        {
            // Define search bracket
            scalar endA =
            (
                (2.0 * alpha) - alphaPrev
            );

            scalar endB =
            (
                Foam::min(alphaMax, alpha + t1 * (alpha - alphaPrev))
            );

            // Find new alpha within bracket
            scalar newAlpha =
            (
                minimizeAlpha
                (
                    endA, endB,
                    alphaPrev, alpha,
                    fPrev, dfPrev,
                    fAlpha, dfAlpha
                )
            );

            alphaPrev = alpha;
            alpha = newAlpha;
        }
        else
        {
            alpha = alphaMax;
        }
    }

    // Bail out if acceptable bracket
    // was not found within specified
    // number of iterations
    if (flag != 2)
    {
        return alpha;
    }

    // Step 2: Find an acceptable point within specified bracket
    iter = 0;

    scalar aPrev, faPrev, dfaPrev;
    scalar bPrev, fbPrev, dfbPrev;

    while (iter < maxIter)
    {
        // Choose a reduced bracket
        scalar endA = a + t2 * (b - a);
        scalar endB = b - t3 * (b - a);

        // Obtain new alpha
        alpha =
        (
            minimizeAlpha
            (
                endA, endB,
                a, b,
                fa, dfa,
                fb, dfb
            )
        );

        // Check for convergence problems
        // due to round-off
        if (Foam::mag((alpha - a) * dfa) <= eps_)
        {
            flag = -2;
            return alpha;
        }

        // Evaluate function and gradient
        fAlpha =
        (
            evaluateFunctional
            (
                data.cellIndex(),
                data.fraction(),
                data.refCentre(),
                (x + (alpha * dir)),
                gradAlpha,
                data.centre(),
                data.distance()
            )
        );

        dfAlpha = (gradAlpha & dir);

        iter++;
        fnEvals++;

        // Update brackets
        aPrev = a; faPrev = fa; dfaPrev = dfa;
        bPrev = b; fbPrev = fb; dfbPrev = dfb;

        if (fAlpha > (fInit + alpha * rho * dfInit) || fAlpha >= fa)
        {
            a = aPrev; b = alpha;
            fa = faPrev; fb = fAlpha;
            dfa = dfaPrev; dfb = dfAlpha;
        }
        else
        {
            // Check if point is acceptable
            if (mag(dfAlpha) <= -sigma * dfInit)
            {
                flag = 0;
                return alpha;
            }

            // Update for next iterate
            a = alpha; fa = fAlpha; dfa = dfAlpha;

            if ((b - a) * dfAlpha >= 0.0)
            {
                b = aPrev; fb = faPrev; dfb = dfaPrev;
            }
            else
            {
                b = bPrev; fb = fbPrev; dfb = dfbPrev;
            }
        }

        // Check for round-off
        if (mag(b - a) < eps_)
        {
            flag = -2;
            return alpha;
        }
    }

    // At this point, the
    // max number of iterations is reached
    flag = -1;
    return alpha;
}


// Broyden-Fletcher-Goldfarb-Shanno (BFGS) algorithm
scalar MomentOfFluid::BFGS
(
    vector2D& x,
    optInfo& data
)
{
    vector2D grad;
    label flag = -1, fnEvals = 0;
    label iter = 0, maxIter = 200;

    // Initial function evaluation
    scalar f =
    (
        evaluateFunctional
        (
            data.cellIndex(),
            data.fraction(),
            data.refCentre(),
            x,
            grad,
            data.centre(),
            data.distance()
        )
    );

    fnEvals++;

    // Infinity norm of initial gradient
    scalar gNorm =
    (
        Foam::max
        (
            Foam::mag(grad[0]),
            Foam::mag(grad[1])
        )
    );

    // Line-search parameters
    scalar rho = 0.01, sigma = 0.9;
    scalar fMin = f - 1e+08 * (1 + mag(f));

    // Specify tolerances
    scalar xTol = 1e-06, fTol = 1e-06;

    scalar fOld, dDir;
    vector2D dx, dir, gradOld;
    scalar alpha = 0.0, initAlpha, alphaOld;

    // Specify initial Hessian
    tensor2D I
    (
        1.0, 0.0,
        0.0, 1.0
    );

    tensor2D H = I;

    // Check for initial convergence
    bool done = (gNorm < fTol);

    while (!done)
    {
        iter++;

        // Specify search direction
        dir = -(H & grad);
        dDir = (grad & dir);

        // Line-search along direction
        initAlpha = 1.0;

        if (iter == 1)
        {
            initAlpha = Foam::min(1.0, (1.0 / gNorm));
        }

        // Store old values
        fOld = f;
        gradOld = grad;
        alphaOld = alpha;

        // Execute line-search
        alpha =
        (
            lineSearch
            (
                x, dir,
                f, dDir,
                initAlpha, rho, sigma,
                fMin, f, grad,
                data, flag, fnEvals
            )
        );

        // If line-search failed, reset and break out
        if (flag < 0 && f >= fOld)
        {
            f = fOld;
            grad = gradOld;
            alpha = alphaOld;

            break;
        }

        // Update using alpha
        dx = (alpha * dir);
        x += dx;

        if (debug)
        {
            Info<< " Iteration: " << iter
                << " fnEvals: " << fnEvals
                << " function: " << f
                << " step-size: " << alpha
                << " norm: " << max(mag(grad[0]), mag(grad[1]))
                << endl;
        }

        // Check if we're done
        scalar gCheck =
        (
            Foam::max
            (
                Foam::mag(grad[0]),
                Foam::mag(grad[1])
            )
        );

        scalar xCheck =
        (
            Foam::max
            (
                Foam::mag(dx[0] / (1.0 + mag(x[0]))),
                Foam::mag(dx[1] / (1.0 + mag(x[1])))
            )
        );

        if
        (
            gCheck < (fTol * (1.0 + gNorm)) ||
            xCheck < xTol ||
            flag == -2 ||
            iter > maxIter
        )
        {
            done = true;
        }
        else
        {
            // Perform BFGS update of Hessian
            vector2D dGrad = (grad - gradOld);
            scalar dXdGrad = (dx & dGrad);

            // Check if update is possible
            bool update =
            (
                dXdGrad >= (sqrteps_ * Foam::max(eps_, mag(dx) * mag(dGrad)))
            );

            if (update)
            {
                if (iter == 1)
                {
                    H = dXdGrad / (dGrad & dGrad) * I;
                }

                vector2D HdGrad = (H & dGrad);

                H +=
                (
                    (1.0 + (dGrad & HdGrad) / dXdGrad) * (dx * dx) / dXdGrad
                  - (dx * HdGrad + HdGrad * dx) / dXdGrad
                );
            }
        }
    }

    if (iter >= maxIter)
    {
        Info<< " Max iterations reached: "
            << "   Iteration: " << iter << nl
            << "   fnEvals: " << fnEvals << nl
            << "   function: " << f << nl
            << "   step-size: " << alpha << nl
            << "   norm: " << max(mag(grad[0]), mag(grad[1])) << nl
            << endl;
    }

    return f;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

MomentOfFluid::MomentOfFluid
(
    const polyMesh& mesh
)
:
    mesh_(mesh),
    tetDecomp_(10),
    allTets_(10),
    allTris_(10)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

MomentOfFluid::~MomentOfFluid()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void MomentOfFluid::constructInterface
(
    const scalarField& fractions,
    const vectorField& refCentres
)
{
    vector normal = vector::zero;
    vector centre = vector::zero;

    scalar minBound = 0.0, maxBound = 1.0;

    forAll(fractions, cellI)
    {
        scalar fraction = fractions[cellI];

        if (fraction > minBound && fraction < maxBound)
        {
            optimizeCentroid
            (
                cellI,
                fraction,
                refCentres[cellI],
                normal,
                centre
            );
        }
    }
}


// Output trianglulated surface to VTK
void MomentOfFluid::outputSurface() const
{
    // Make the directory
    fileName dirName(mesh_.time().path()/"VTK"/mesh_.time().timeName());

    mkDir(dirName);

    // Open stream for output
    OFstream file(dirName/"MoF.vtk");

    label nTris = allTris_.size();
    label nPoints = (3 * nTris);

    // Write out the header
    file<< "# vtk DataFile Version 2.0" << nl
        << "MoF.vtk" << nl
        << "ASCII" << nl
        << "DATASET UNSTRUCTURED_GRID" << nl
        << "POINTS " << nPoints << " double" << nl;

    forAll(allTris_, triI)
    {
        const MoF::Triangle& tri = allTris_[triI];

        forAll(tri, i)
        {
            file<< tri[i].x() << ' '
                << tri[i].y() << ' '
                << tri[i].z() << ' '
                << nl;
        }
    }

    file << "CELLS " << nTris << " " << (nPoints + nTris) << endl;

    nPoints = 0;

    forAll(allTris_, triI)
    {
        const MoF::Triangle& tri = allTris_[triI];

        file<< 3 << ' ';

        forAll(tri, i)
        {
            file<< nPoints++ << ' ';
        }

        file<< nl;
    }

    file<< "CELL_TYPES " << nTris << endl;

    forAll(allTris_, triI)
    {
        file<< 5 << ' ';
    }
}


// Output plane as VTK
void MomentOfFluid::outputPlane
(
    const point& p,
    const vector& n,
    const label& cellIndex
) const
{
    // Make the directory
    fileName dirName(mesh_.time().path()/"VTK"/mesh_.time().timeName());

    mkDir(dirName);

    // Open stream for output
    OFstream file(dirName/(Foam::name(n)+".vtk"));

    FixedList<point, 4> planePoints;

    // Configure points on plane
    vector norm = n / (mag(n) + VSMALL);

    const cell& checkCell = mesh_.cells()[cellIndex];

    pointField pts = checkCell.points(mesh_.faces(), mesh_.points());

    scalar maxMag = 0.0;
    vector R = vector::zero;
    vector T = vector::zero;

    forAll(pts, pI)
    {
        vector r = pts[pI] - p;

        if (mag(r) > maxMag)
        {
            R = r;
            maxMag = mag(r);
        }
    }

    // Scale it up a little bit
    R *= 1.5;

    // Remove normal component
    R = R - ((R & norm) * norm);
    T = (norm ^ R);

    planePoints[0] = p + R;
    planePoints[1] = p + T;
    planePoints[2] = p - R;
    planePoints[3] = p - T;

    // Write out the header
    file<< "# vtk DataFile Version 2.0" << nl
        << "MoF.vtk" << nl
        << "ASCII" << nl
        << "DATASET UNSTRUCTURED_GRID" << nl
        << "POINTS 4 double" << nl;

    forAll(planePoints, pI)
    {
        file<< planePoints[pI].x() << ' '
            << planePoints[pI].y() << ' '
            << planePoints[pI].z() << ' '
            << nl;
    }

    file<< "CELLS 1 5" << endl;
    file<< "4 0 1 2 3" << endl;
    file<< "CELL_TYPES 1" << endl;
    file<< "9" << endl;
}


} // End namespace Foam

// ************************************************************************* //
