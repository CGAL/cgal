// Copyright (c) 2006-2009 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Kerber <mkerber@mpi-inf.mpg.de>
//
// ============================================================================



#ifndef CGAL_ACK_FLAGS_H
#define CGAL_ACK_FLAGS_H 1

// TODO: use new interface and remove this flag
#ifndef CGAL_AK_ENABLE_DEPRECATED_INTERFACE
#define CGAL_AK_ENABLE_DEPRECATED_INTERFACE 0
#endif

// Is debug-information printed?
#ifndef CGAL_ACK_DEBUG_FLAG
#define CGAL_ACK_DEBUG_FLAG 0
#endif

// If CGAL_ACK_DEBUG_FLAG is set, which output stream is used for debug?
#ifndef CGAL_ACK_DEBUG_PRINT
#define CGAL_ACK_DEBUG_PRINT std::cout
#endif

// If enabled, needs includes from experimental package
#ifndef CGAL_ACK_WITH_FILTERED_KERNEL
#define CGAL_ACK_WITH_FILTERED_KERNEL 0
#endif

// If enabled, needs includes from experimental package
#ifndef CGAL_ACK_WITH_ROTATIONS
#define CGAL_ACK_WITH_ROTATIONS 0
#endif

/**
 * The threshold that is used in the Filtered_algebraic_curve_kernel_2
 */
#ifndef CGAL_ACK_THRESHOLD_FOR_FILTERED_KERNEL
#define CGAL_ACK_THRESHOLD_FOR_FILTERED_KERNEL 0.01
#endif

/**
 * For random choices in the algorithm, this seed is used
 * If set to zero, a random seed is used
 */
#ifndef CGAL_ACK_STATIC_SEED
#define CGAL_ACK_STATIC_SEED 0
#endif

/**
 * Allows to use the Bitstream tree described in Eigenwillig's thesis
 */
#ifndef CGAL_ACK_BITSTREAM_USES_E08_TREE
#define CGAL_ACK_BITSTREAM_USES_E08_TREE 1
#endif


/**
 * If set, the program uses the AlciX-code
 * for the curve- and curve-pair-analysis.
 * This flag is only for debugging purposes.
 */
#ifndef CGAL_ACK_USE_EXACUS
#define CGAL_ACK_USE_EXACUS 0
#endif

/**
 * If set, the curve and curve pair analysis are using specialized code
 * to analyse conic curves, i.e. curves of degree 2
 */
#ifndef CGAL_ACK_USE_SPECIAL_TREATMENT_FOR_CONIX
#define CGAL_ACK_USE_SPECIAL_TREATMENT_FOR_CONIX 0
#endif

/**
 * The curve analysis does not distinguish between "(1,1)-singularities"
 * (i.e., vertical cusps, isolated points on arcs), and usual regular points.
 * The candidate point on each status line can be checked for being singular
 * using this flag. This gives additional information but increases
 * compuation time
 *
 * WARNING: Currently, the status line does not store the additional
 * information whether a point is singluar or not.
 * Therefore, there is currently no reasons to set this flag. It is still
 * contained for possible further extension of the status line.
 */
#ifndef CGAL_ACK_CHECK_CANDIDATE_FOR_SINGULARITY
#define CGAL_ACK_CHECK_CANDIDATE_FOR_SINGULARITY 0
#endif

/**
 * If set to 1, curve pairs are not checked for coprimality. Only do this
 * if you know what you are doing!
 */
#ifndef CGAL_ACK_DONT_CHECK_POLYNOMIALS_FOR_COPRIMALITY
#define CGAL_ACK_DONT_CHECK_POLYNOMIALS_FOR_COPRIMALITY 0
#endif

/**
 * The "resultant first" strategy means: instead of computing the full
 * subresultant sequence (or Sturm-Habicht sequence), the algorithm
 * only computed the resultant in a first step. This suffices already for
 * many curves (i.e., regular ones). The full subresultant is computed
 * if it is needed for the first time.
 *
 * This strategy only makes sense if computing resultants is faster than
 * computing subresultants, otherwise, it wastes computation time.
 * Since resultant computation is done by interpolation,
 * it is faster than the pseudo-remainder based subresultant computation.
 */
#ifndef CGAL_ACK_RESULTANT_FIRST_STRATEGY
#define CGAL_ACK_RESULTANT_FIRST_STRATEGY 1
#endif

/**
 * If CGAL_ACK_RESULTANT_FIRST_STRATEGY is set, this flag determines
 * for which curves the "resultant first" strategy is used. Depending
 * on the resultant algorithm, the strategy might only be advantageous
 * for higher degree curves
 *
 * If CGAL_ACK_RESULTANT_FIRST_STRATEGY is not set, this flag has no effect
 */
#ifndef CGAL_ACK_RESULTANT_FIRST_STRATEGY_DEGREE_THRESHOLD
#define CGAL_ACK_RESULTANT_FIRST_STRATEGY_DEGREE_THRESHOLD 0
#endif

/**
 * Subresultants can be computed by a polynomial remainder sequence (default),
 * or by evaluating minors of the bezout matrix by setting this flag.
 * Tests have shown that the polynomial remainder sequence is more efficient
 * so it is not recommended to set this flag.
 */
#ifndef CGAL_ACK_USE_BEZOUT_MATRIX_FOR_SUBRESULTANTS
#define CGAL_ACK_USE_BEZOUT_MATRIX_FOR_SUBRESULTANTS 0
#endif

/**
 * Allows to switch off the specialized method for Status_line_CPA_1
 * if multiplicity is zero or one.
 * Since this methods improves the performance,
 * it is not recommended to set this flag unless for testing
 */
#ifndef ACK_CGAL_NO_ARC_FLIP
#define ACK_CGAL_NO_ARC_FLIP 0
#endif

/**
 * This flags defines the default strategy to handle degenerate curves
 * There are three choices currently available:
 * SHEAR_STRATEGY performs a shear whenever a degenerate situation occurs.
 * SHEAR_ONLY_AT_IRRATIONAL_STRATEGY handles rational coordinates with
 * a more direct method, but performs a shear for irrational x-coordinates
 * that have a degeneracy. Finally, EXCEPTION_STRATEGY throws an exception
 * whenever a degeneracy occurs.
 */
#ifndef CGAL_ACK_DEFAULT_DEGENERACY_STRATEGY
//#define CGAL_ACK_DEFAULT_DEGENERACY_STRATEGY CGAL::SHEAR_ONLY_AT_IRRATIONAL_STRATEGY
#define CGAL_ACK_DEFAULT_DEGENERACY_STRATEGY CGAL::SHEAR_STRATEGY
#endif

/**
 * The algorithm can also handle non-y-regular curves without shearing,
 * in case that the resultant multiplicity at vertical asymptotes is one.
 * This special treatement can be switched off by setting this flag.
 * It is not recommended to do this because of efficiency
 */
#ifndef CGAL_ACK_SHEAR_ALL_NOT_Y_REGULAR_CURVES
#define CGAL_ACK_SHEAR_ALL_NOT_Y_REGULAR_CURVES 0
#endif

/**
 * At some points in the algorithm, it is checked whether a polynomial
 * H(x):=h(p(x),q(x)) vanishes for an algebraic number x_0 with polynomial r.
 * For that check, the computation of H is done modulo r for efficiency.
 * This can be switched off by this flag, though it is recommended not to
 * do so.
 */
#ifndef CGAL_ACK_USE_NO_REDUCTION_MODULO_RESULTANT
#define CGAL_ACK_USE_NO_REDUCTION_MODULO_RESULTANT 0
#endif

#ifndef CGAL_AK_DONT_USE_SIMPLE_BOUND_BETWEEN
#define CGAL_AK_DONT_USE_SIMPLE_BOUND_BETWEEN 0
#endif

/**
 * These flags are experimental, concerning interval arithmetic methods
 * Don't change them!
 */
#ifndef CGAL_ACK_USE_DERIVATIVE_OPTION
#define CGAL_ACK_USE_DERIVATIVE_OPTION 0
#endif
#ifndef CGAL_ACK_USE_BISECTION_OPTION
#define CGAL_ACK_USE_BISECTION_OPTION 0
#endif


#endif // CGAL_ACK_FLAGS_H
