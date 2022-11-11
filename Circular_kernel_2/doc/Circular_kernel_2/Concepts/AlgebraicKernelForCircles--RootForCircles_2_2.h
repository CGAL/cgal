
/*!
\ingroup PkgCircularKernel2AlgebraicConcepts
\cgalConcept

Concept to represent the roots of a system of two equations of degree 2
in two variables `x` and `y` that are models of concept
`AlgebraicKernelForCircles::PolynomialForCircles_2_2`

\cgalHasModel CGAL::Root_for_circles_2_2

\sa `AlgebraicKernelForCircles`

*/

class AlgebraicKernelForCircles::RootForCircles_2_2 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!

*/
bool operator ==( const AlgebraicKernelForCircles::RootForCircles_2_2 & p,
                  const AlgebraicKernelForCircles::RootForCircles_2_2 & q);

/// @}

}; /* end AlgebraicKernelForCircles::RootForCircles_2_2 */

