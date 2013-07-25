
namespace CGAL {

/*!
\ingroup PkgCircularKernel2GeometricClasses

\cgalModels `CircularKernel`

\cgalHeading{Parameters}

The first parameter of the circular kernel must be instantiated with a 
model of the `Kernel` concept. The `Circular_kernel_2` class 
template derives from this first parameter, in order to reuse all 
needed functionalities on basic objects provided by a model of the 
Kernel concept. 

The second parameter, `AlgebraicKernelForCircles`, is meant to provide the 
circular kernel with all the algebraic functionalities required for the 
manipulation of algebraic curves. 

\sa `Kernel`
\sa `AlgebraicKernelForCircles`
\sa `CGAL::Exact_circular_kernel_2`

*/
template< typename Kernel, typename AlgebraicKernelForCircles >
class Circular_kernel_2 : public Kernel {
public:

/// \name Types 
/// The circular kernel uses basic number types of the algebraic
/// kernel. In fact, the two number types
/// `AlgebraicKernelForCircles::RT` and `Kernel::RT` must coincide, as
/// well as `AlgebraicKernelForCircles::FT` and `Kernel::FT`.
/// @{

/*!
Ring number type. 
*/ 
typedef AlgebraicKernelForCircles::RT RT; 

/*!
Field number type. 
*/ 
typedef AlgebraicKernelForCircles::FT FT; 

/// @}

/// \name
/// The following types are available, as well as all the
/// functionality on them described in the `CircularKernel` concept.
/// @{

/*!

*/ 
typedef Line_arc_2<Circular_kernel_2> Line_arc_2; 

/*!

*/ 
typedef Circular_arc_2<Circular_kernel_2> Circular_arc_2; 

/*!

*/ 
typedef Circular_arc_point_2<Circular_kernel_2> Circular_arc_point_2; 

/// @}

}; /* end Circular_kernel_2 */
} /* end namespace CGAL */
