
/*!
\ingroup PkgSphericalKernel3GeometricConcepts
\cgalconcept

An object `fo` of this type must provide: 

*/

class SphericalKernel::ComputeApproximateAngle_3 {
public:

/*! 
Computes an approximation of the angle of the arc in radians `a`. 
*/ 
double 
operator()(const SphericalKernel::Circular_arc_3 & a); 
}; /* end SphericalKernel::ComputeApproximateAngle_3 */

