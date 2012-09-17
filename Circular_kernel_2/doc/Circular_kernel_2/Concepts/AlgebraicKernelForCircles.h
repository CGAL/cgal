
/*!
\ingroup PkgCircularKernel2Concepts
\cgalconcept

The `AlgebraicKernelForCircles` concept is meant to provide the 
curved kernel with all the algebraic functionalities required for the 
manipulation of circular arcs. 

\hasModel Algebraic_kernel_for_circles_2_2 

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
typedef Hidden_type RT; 

/*! 
A model of `FieldNumberType``<RT>`. 
*/ 
typedef Hidden_type FT; 

/*! 
A model of 
`AlgebraicKernelForCircles::Polynomial_1_2`, for bivariate polynomials of degree up 
to 1. 
*/ 
typedef Hidden_type Polynomial_1_2; 

/*! 
A model of 
`AlgebraicKernelForCircles::PolynomialForCircles_2_2`, for bivariate polynomials 
of degree up to 2 that can store equations of circles. 
*/ 
typedef Hidden_type Polynomial_for_circles_2_2; 

/*! 
A model of 
`RootOf_2`, for algebraic numbers 
of degree up to 2. 
*/ 
typedef Hidden_type Root_of_2; 

/*! 
A model of 
`AlgebraicKernelForCircles::RootForCircles_2_2`, for 
solutions of systems of two models of 
`AlgebraicKernelForCircles::PolynomialForCircles_2_2`. 
*/ 
typedef Hidden_type Root_for_circles_2_2; 

/*! 
A model of 
`AlgebraicKernelForCircles::ConstructPolynomial_1_2`. 
*/ 
typedef Hidden_type Construct_polynomial_1_2; 

/*! 
A model of 
`AlgebraicKernelForCircles::ConstructPolynomialForCircles_2_2`. 
*/ 
typedef Hidden_type Construct_polynomial_for_circles_2_2; 

/*! 
A model of the concept 
`AlgebraicKernelForCircles::CompareX`. 
*/ 
typedef Hidden_type Compare_x; 

/*! 
A model of the concept 
`AlgebraicKernelForCircles::CompareY`. 
*/ 
typedef Hidden_type Compare_y; 

/*! 
A model of the concept 
`AlgebraicKernelForCircles::CompareXY`. 
*/ 
typedef Hidden_type Compare_xy; 

/*! 
A model of the concept `AlgebraicKernelForCircles::SignAt`. 
*/ 
typedef Hidden_type Sign_at; 

/*! 
A model of the concept 
`AlgebraicKernelForCircles::XCriticalPoints`. 
*/ 
typedef Hidden_type X_critical_points; 

/*! 
A model of the concept 
`AlgebraicKernelForCircles::YCriticalPoints`. 
*/ 
typedef Hidden_type Y_critical_points; 

/*! 
A model of the concept `AlgebraicKernelForCircles::Solve`. 
*/ 
typedef Hidden_type Solve; 

/// @}

}; /* end AlgebraicKernelForCircles */

