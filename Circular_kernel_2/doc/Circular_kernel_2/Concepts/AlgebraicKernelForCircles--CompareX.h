
/*!
\ingroup PkgCircularKernel2AlgebraicConcepts
\cgalConcept

\sa `AlgebraicKernelForCircles::CompareY`
\sa `AlgebraicKernelForCircles::CompareXY`
\sa `CircularKernel::CompareX_2`

*/

class AlgebraicKernelForCircles::CompareX {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
Compares the `x` (first) variables of two `Root_for_circles_2_2`.
*/
template < class OutputIterator >
CGAL::Comparison_result
operator()(const AlgebraicKernelForCircles::Root_for_circles_2_2 & r1,
const AlgebraicKernelForCircles::Root_for_circles_2_2 & r2);

/// @}

}; /* end AlgebraicKernelForCircles::CompareX */

