
/*!
\ingroup PkgSphericalKernel3Concepts
\cgalconcept

The `AlgebraicKernelForSpheres` concept is meant to provide the 
curved kernel with all the algebraic functionalities required for the 
manipulation of spheres, circles, and circular arcs in 3D. 

\hasModel Algebraic_kernel_for_spheres_2_3 

\sa `SphericalKernel`
\sa `CGAL::Spherical_kernel_3<Kernel,AlgebraicKernelForSpheres>`

*/

class AlgebraicKernelForSpheres {
public:

/// \name Types 
/// A model of `AlgebraicKernelForSpheres` is supposed to provide
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
`AlgebraicKernelForSpheres::Polynomial_1_3`, for trivariate polynomials 
of degree up to 1. 
*/ 
typedef Hidden_type Polynomial_1_3; 

/*! 
A model of 
`AlgebraicKernelForSpheres::PolynomialForSpheres_2_3`, for trivariate 
polynomials of degree up to 2 that can store equations of spheres. 
*/ 
typedef Hidden_type Polynomial_for_spheres_2_3; 

/*! 
A model of 
`AlgebraicKernelForSpheres::PolynomialsForLines_3`, for systems of 
polynomials that can store equations of lines in 3D. 
*/ 
typedef Hidden_type Polynomials_for_lines_3; 

/*! 
A model of 
`RootOf_2`, for algebraic numbers 
of degree up to 2. 
*/ 
typedef Hidden_type Root_of_2; 

/*! 
A model of 
`AlgebraicKernelForSpheres::RootForSpheres_2_3`, for 
solutions of systems of three models of 
`AlgebraicKernelForSpheres::PolynomialForSpheres_2_3`. 
*/ 
typedef Hidden_type Root_for_spheres_2_3; 

/*! 
A model of 
`AlgebraicKernelForSpheres::ConstructPolynomial_1_3`. 
*/ 
typedef Hidden_type Construct_polynomial_1_3; 

/*! 
A model of 
`AlgebraicKernelForSpheres::ConstructPolynomialForSpheres_2_3`. 
*/ 
typedef Hidden_type Construct_polynomial_for_spheres_2_3; 

/*! 
A model of 
`AlgebraicKernelForSpheres::ConstructPolynomialsForLines_3`. 
*/ 
typedef Hidden_type Construct_polynomials_for_lines_3; 

/*! 
A model of the concept 
`AlgebraicKernelForSpheres::CompareX`. 
*/ 
typedef Hidden_type Compare_x; 

/*! 
A model of the concept 
`AlgebraicKernelForSpheres::CompareY`. 
*/ 
typedef Hidden_type Compare_y; 

/*! 
A model of the concept 
`AlgebraicKernelForSpheres::CompareZ`. 
*/ 
typedef Hidden_type Compare_z; 

/*! 
A model of the concept 
`AlgebraicKernelForSpheres::CompareXY`. 
*/ 
typedef Hidden_type Compare_xy; 

/*! 
A model of the concept 
`AlgebraicKernelForSpheres::CompareXYZ`. 
*/ 
typedef Hidden_type Compare_xyz; 

/*! 
A model of the concept `AlgebraicKernelForSpheres::SignAt`. 
*/ 
typedef Hidden_type Sign_at; 

/*! 
A model of the concept 
`AlgebraicKernelForSpheres::XCriticalPoints`. 
*/ 
typedef Hidden_type X_critical_points; 

/*! 
A model of the concept 
`AlgebraicKernelForSpheres::YCriticalPoints`. 
*/ 
typedef Hidden_type Y_critical_points; 

/*! 
A model of the concept 
`AlgebraicKernelForSpheres::ZCriticalPoints`. 
*/ 
typedef Hidden_type Z_critical_points; 

/*! 
A model of the concept `AlgebraicKernelForSpheres::Solve`. 
*/ 
typedef Hidden_type Solve; 

/// @}

}; /* end AlgebraicKernelForSpheres */

