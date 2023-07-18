
/*!
\ingroup PkgCircularKernel2AlgebraicConcepts
\cgalConcept

The `AlgebraicKernelForCircles` concept is meant to provide the
curved kernel with all the algebraic functionalities required for the
manipulation of circular arcs.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Algebraic_kernel_for_circles_2_2}
\cgalHasModelsEnd

\sa `CircularKernel`
\sa `CGAL::Circular_kernel_2<Kernel,AlgebraicKernelForCircles>`

*/

class AlgebraicKernelForCircles {
public:

/// \name Types
/// A model of `AlgebraicKernelForCircles` is supposed to provide
/// @{

/*!
A model of `RingNumberType`.
*/
typedef unspecified_type RT;

/*!
A model of `FieldNumberType``<RT>`.
*/
typedef unspecified_type FT;

/*!
A model of
`AlgebraicKernelForCircles::Polynomial_1_2`, for bivariate polynomials of degree up
to 1.
*/
typedef unspecified_type Polynomial_1_2;

/*!
A model of
`AlgebraicKernelForCircles::PolynomialForCircles_2_2`, for bivariate polynomials
of degree up to 2 that can store equations of circles.
*/
typedef unspecified_type Polynomial_for_circles_2_2;

/*!
A model of
`RootOf_2`, for algebraic numbers
of degree up to 2.
*/
typedef unspecified_type Root_of_2;

/*!
A model of
`AlgebraicKernelForCircles::RootForCircles_2_2`, for
solutions of systems of two models of
`AlgebraicKernelForCircles::PolynomialForCircles_2_2`.
*/
typedef unspecified_type Root_for_circles_2_2;

/*!
A model of
`AlgebraicKernelForCircles::ConstructPolynomial_1_2`.
*/
typedef unspecified_type Construct_polynomial_1_2;

/*!
A model of
`AlgebraicKernelForCircles::ConstructPolynomialForCircles_2_2`.
*/
typedef unspecified_type Construct_polynomial_for_circles_2_2;

/*!
A model of the concept
`AlgebraicKernelForCircles::CompareX`.
*/
typedef unspecified_type Compare_x;

/*!
A model of the concept
`AlgebraicKernelForCircles::CompareY`.
*/
typedef unspecified_type Compare_y;

/*!
A model of the concept
`AlgebraicKernelForCircles::CompareXY`.
*/
typedef unspecified_type Compare_xy;

/*!
A model of the concept `AlgebraicKernelForCircles::SignAt`.
*/
typedef unspecified_type Sign_at;

/*!
A model of the concept
`AlgebraicKernelForCircles::XCriticalPoints`.
*/
typedef unspecified_type X_critical_points;

/*!
A model of the concept
`AlgebraicKernelForCircles::YCriticalPoints`.
*/
typedef unspecified_type Y_critical_points;

/*!
A model of the concept `AlgebraicKernelForCircles::Solve`.
*/
typedef unspecified_type Solve;

/// @}

}; /* end AlgebraicKernelForCircles */

