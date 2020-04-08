
/*!
\ingroup PkgCircularKernel2AlgebraicConcepts
\cgalConcept

\sa `AlgebraicKernelForCircles::YCriticalPoints`

*/

class AlgebraicKernelForCircles::XCriticalPoints {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
Copies in the output iterator the `x`-critical points of polynomial
`p`, as objects of type `AlgebraicKernelForCircles::Root_for_circles_2_2`.
*/
template < class OutputIterator >
OutputIterator
operator()(const AlgebraicKernelForCircles::Polynomial_for_circles_2_2 &p,
OutputIterator res);

/*!
Computes the `x`-critical point with smallest (resp.\ largest) `x`
of polynomial `p` if `b` is `true` (resp.\ `false`).
*/
template < class OutputIterator >
AlgebraicKernelForCircles::Root_for_circles_2_2
operator()(const AlgebraicKernelForCircles::Polynomial_for_circles_2_2 &p,
bool b);

/// @}

}; /* end AlgebraicKernelForCircles::XCriticalPoints */

