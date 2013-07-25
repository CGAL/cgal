
/*!
\ingroup PkgCircularKernel2AlgebraicConcepts
\cgalConcept

\sa `AlgebraicKernelForCircles::XCriticalPoints`

*/

class AlgebraicKernelForCircles::YCriticalPoints {
public:

/// \name Operations
/// A model `fo` of this type must provide:
/// @{

/*!
Copies in the output iterator the `y`-critical points of polynomial 
`p`, as objects of type `AlgebraicKernelForCircles::Root_for_circles_2_2`. 
*/ 
template < class OutputIterator > 
OutputIterator 
operator()(const AlgebraicKernelForCircles::Polynomial_for_circles_2_2 &p, 
OutputIterator res); 

/*!
Computes the `y`-critical point with smallest (resp.\ largest) `y` 
of polynomial `p` if `b` is `true` (resp.\ `false`). 
*/ 
template < class OutputIterator > 
AlgebraicKernelForCircles::Root_for_circles_2_2 
operator()(const AlgebraicKernelForCircles::Polynomial_for_circles_2_2 &p, 
bool i); 

/// @}

}; /* end AlgebraicKernelForCircles::YCriticalPoints */

