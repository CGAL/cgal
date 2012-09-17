
/*!
\ingroup PkgCircularKernel2Concepts
\cgalconcept

A model `fo` of this type must provide: 

\sa `AlgebraicKernelForCircles::CompareY`
\sa `AlgebraicKernelForCircles::CompareXY`
\sa `CircularKernel::CompareX_2`

*/

class AlgebraicKernelForCircles::CompareX {
public:

/// \name See Also 
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

