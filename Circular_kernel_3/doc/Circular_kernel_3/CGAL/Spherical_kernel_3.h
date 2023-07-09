
namespace CGAL {

/*!
\ingroup PkgCircularKernel3GeometricClasses

\cgalModels{SphericalKernel}

\cgalHeading{Parameters}

The first parameter of the spherical kernel must be instantiated with
a model of the `Kernel` concept. The `Spherical_kernel_3`
class template derives from this first parameter, in order to reuse
all needed functionalities on basic objects provided by a model of the
Kernel concept.

The second parameter, `AlgebraicKernelForSpheres`, is meant to provide the
spherical kernel with all the algebraic functionalities required for the
manipulation of algebraic curves.

\sa `Kernel`
\sa `AlgebraicKernelForSpheres`
\sa `CGAL::Exact_spherical_kernel_3`

*/
template< typename Kernel, typename AlgebraicKernelForSpheres >
struct Spherical_kernel_3 : public Kernel {

/// \name Types
/// The spherical kernel uses basic number types of the algebraic
/// kernel: In fact, the two number types
/// `AlgebraicKernelForSpheres::RT` and `Kernel::RT` must refer to the
/// same type, as well as `AlgebraicKernelForSpheres::FT` and
/// `Kernel::FT`.
/// @{

/*!
Ring number type.
*/
typedef AlgebraicKernelForSpheres::RT RT;

/*!
Field number type.
*/
typedef AlgebraicKernelForSpheres::FT FT;

/// @}

/// \name
/// The following types are available, as well as all the
/// functionality on them described in the `SphericalKernel`
/// concept.
///
/// `Polynomials_for_circle_3` is implemented as a
/// `std::pair< Polynomial_for_spheres_2_3, Polynomial_1_3 >` and is a
/// model of `AlgebraicKernelForSpheres::PolynomialsForCircles_3`.
/// @{

/*!

*/
typedef Line_arc_3<Spherical_kernel_3> Line_arc_3;

/*!

*/
typedef Circular_arc_3<Spherical_kernel_3> Circular_arc_3;

/*!

*/
typedef Circular_arc_point_3<Spherical_kernel_3> Circular_arc_point_3;

/// @}

}; /* end Spherical_kernel_3 */
} /* end namespace CGAL */
