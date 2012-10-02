
/*!
\ingroup PkgSphericalKernel3GeometricConcepts
\cgalconcept

\refines ::Kernel::ComputeApproximateSquaredLength_3 
In addition, an object `fo` of this type must provide: 

*/

class SphericalKernel::ComputeApproximateSquaredLength_3 {
public:
/*! 
Computes an approximation of the squared length of the arc `a`. 
*/ 
double 
operator()(const SphericalKernel::Circular_arc_3 & a); 

}; /* end SphericalKernel::ComputeApproximateSquaredLength_3 */

