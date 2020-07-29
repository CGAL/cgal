
/*!
\ingroup PkgCircularKernel2AlgebraicConcepts
\cgalConcept

\sa `AlgebraicKernelForCircles::CompareX`
\sa `AlgebraicKernelForCircles::CompareXY`
\sa `CircularKernel::CompareY_2`

*/

class AlgebraicKernelForCircles::CompareY {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
Compares the `y` (second) variables of two `Root_for_circles_2_2`.
*/
template < class OutputIterator >
CGAL::Comparison_result
operator()(const AlgebraicKernelForCircles::Root_for_circles_2_2 & r1,
const AlgebraicKernelForCircles::Root_for_circles_2_2 & r2);

/// @}

}; /* end AlgebraicKernelForCircles::CompareY */

