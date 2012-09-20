
/*!
\ingroup PkgSphericalKernel3AlgebraicConcepts
\cgalconcept

Concept to represent the roots of a system of three equations of degree 2 
in three variables `x`, `y` and `z` that are models of concept 
`AlgebraicKernelForSpheres::PolynomialForSpheres_2_3`. 

\hasModel CGAL::Root_for_spheres_2_3 

\sa `AlgebraicKernelForSpheres`

*/

class AlgebraicKernelForSpheres::RootForSpheres_2_3 {
public:

/*! 
Test equality of two roots of systems. 
*/ 
bool operator ==(const AlgebraicKernelForSpheres::RootForSpheres_2_3 & r1, 
                 const AlgebraicKernelForSpheres::RootForSpheres_2_3 & r2); 

}; /* end AlgebraicKernelForSpheres::RootForSpheres_2_3 */

