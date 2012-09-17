
/*!
\ingroup PkgCircularKernel2AlgebraicConcepts
\cgalconcept

Concept to represent 
bivariate polynomials of degree up to 2 capable of storing equations 
of circles, whose center's coordinates, as well as the square of the radius, 
are of a type that is a model of the concept 
`FieldNumberType`. 

\refines CopyConstructible
\refines Assignable
\refines DefaultConstructible 

\hasModel CGAL::Polynomial_for_circles_2_2 

\sa `AlgebraicKernelForCircles`

*/

class AlgebraicKernelForCircles::PolynomialForCircles_2_2 {
public:

/*! 

*/ 
bool operator ==(AlgebraicKernelForCircles::const PolynomialForCircles_2_2 & p, 
                 const AlgebraicKernelForCircles::PolynomialForCircles_2_2 & q); 
}; /* end AlgebraicKernelForCircles::PolynomialForCircles_2_2 */

