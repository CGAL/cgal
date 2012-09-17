
/*!
\ingroup PkgCircularKernel2Concepts
\cgalconcept

Concept to represent the roots of a system of two equations of degree 2 
in two variables `x` and `y` that are models of concept 
`AlgebraicKernelForCircles::PolynomialForCircles_2_2` 

\hasModel CGAL::Root_for_circles_2_2 

\sa `AlgebraicKernelForCircles`

*/

class AlgebraicKernelForCircles::RootForCircles_2_2 {
public:

/*! 

*/ 
bool operator ==( const AlgebraicKernelForCircles::RootForCircles_2_2 & p, 
                  const AlgebraicKernelForCircles::RootForCircles_2_2 & q); 
}; /* end AlgebraicKernelForCircles::RootForCircles_2_2 */

