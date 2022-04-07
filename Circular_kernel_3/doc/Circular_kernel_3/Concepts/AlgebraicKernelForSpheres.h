
/*!
\ingroup PkgCircularKernel3AlgebraicConcepts
\cgalConcept

The `AlgebraicKernelForSpheres` concept is meant to provide the
curved kernel with all the algebraic functionalities required for the
manipulation of spheres, circles, and circular arcs in 3D.

\cgalHasModel CGAL::Algebraic_kernel_for_spheres_2_3

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
typedef unspecified_type RT;

/*!
A model of `FieldNumberType``<RT>`.
*/
typedef unspecified_type FT;

/*!
A model of
`AlgebraicKernelForSpheres::Polynomial_1_3`, for trivariate polynomials
of degree up to 1.
*/
typedef unspecified_type Polynomial_1_3;

/*!
A model of
`AlgebraicKernelForSpheres::PolynomialForSpheres_2_3`, for trivariate
polynomials of degree up to 2 that can store equations of spheres.
*/
typedef unspecified_type Polynomial_for_spheres_2_3;

/*!
A model of
`AlgebraicKernelForSpheres::PolynomialsForLines_3`, for systems of
polynomials that can store equations of lines in 3D.
*/
typedef unspecified_type Polynomials_for_lines_3;

/*!
A model of
`RootOf_2`, for algebraic numbers
of degree up to 2.
*/
typedef unspecified_type Root_of_2;

/*!
A model of
`AlgebraicKernelForSpheres::RootForSpheres_2_3`, for
solutions of systems of three models of
`AlgebraicKernelForSpheres::PolynomialForSpheres_2_3`.
*/
typedef unspecified_type Root_for_spheres_2_3;

/*!
A model of
`AlgebraicKernelForSpheres::ConstructPolynomial_1_3`.
*/
typedef unspecified_type Construct_polynomial_1_3;

/*!
A model of
`AlgebraicKernelForSpheres::ConstructPolynomialForSpheres_2_3`.
*/
typedef unspecified_type Construct_polynomial_for_spheres_2_3;

/*!
A model of
`AlgebraicKernelForSpheres::ConstructPolynomialsForLines_3`.
*/
typedef unspecified_type Construct_polynomials_for_lines_3;

/*!
A model of the concept
`AlgebraicKernelForSpheres::CompareX`.
*/
typedef unspecified_type Compare_x;

/*!
A model of the concept
`AlgebraicKernelForSpheres::CompareY`.
*/
typedef unspecified_type Compare_y;

/*!
A model of the concept
`AlgebraicKernelForSpheres::CompareZ`.
*/
typedef unspecified_type Compare_z;

/*!
A model of the concept
`AlgebraicKernelForSpheres::CompareXY`.
*/
typedef unspecified_type Compare_xy;

/*!
A model of the concept
`AlgebraicKernelForSpheres::CompareXYZ`.
*/
typedef unspecified_type Compare_xyz;

/*!
A model of the concept `AlgebraicKernelForSpheres::SignAt`.
*/
typedef unspecified_type Sign_at;

/*!
A model of the concept
`AlgebraicKernelForSpheres::XCriticalPoints`.
*/
typedef unspecified_type X_critical_points;

/*!
A model of the concept
`AlgebraicKernelForSpheres::YCriticalPoints`.
*/
typedef unspecified_type Y_critical_points;

/*!
A model of the concept
`AlgebraicKernelForSpheres::ZCriticalPoints`.
*/
typedef unspecified_type Z_critical_points;

/*!
A model of the concept `AlgebraicKernelForSpheres::Solve`.
*/
typedef unspecified_type Solve;

/// @}

}; /* end AlgebraicKernelForSpheres */

