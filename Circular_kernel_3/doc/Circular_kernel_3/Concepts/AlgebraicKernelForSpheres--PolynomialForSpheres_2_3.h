
/*!
\ingroup PkgSphericalKernel3Concepts
\cgalconcept

Concept to represent 
trivariate polynomials of degree up to 2 capable of storing equations 
of spheres, whose center's coordinates, as well as the square of the radius, 
are of a type that is a model of the concept 
`FieldNumberType`. 

\refines `CopyConstructible`
\refines `Assignable`
\refines `DefaultConstructible`

\hasModel CGAL::Polynomial_for_spheres_2_3 

\sa `AlgebraicKernelForSpheres`

*/

class AlgebraicKernelForSpheres::PolynomialForSpheres_2_3 {
public:

/*! 
Tests equality of two polynomials. 
*/ 
bool operator ==(const AlgebraicKernelForSpheres::PolynomialForSpheres_2_3 & p, 
const AlgebraicKernelForSpheres::PolynomialForSpheres_2_3 & q); 

}; /* end AlgebraicKernelForSpheres::PolynomialForSpheres_2_3 */

