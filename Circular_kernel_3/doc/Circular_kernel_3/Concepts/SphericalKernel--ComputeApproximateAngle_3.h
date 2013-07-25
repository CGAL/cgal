
/*!
\ingroup PkgSphericalKernel3GeometricConcepts
\cgalConcept
*/

class SphericalKernel::ComputeApproximateAngle_3 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
Computes an approximation of the angle of the arc in radians `a`. 
*/ 
double 
operator()(const SphericalKernel::Circular_arc_3 & a); 

///@}

}; /* end SphericalKernel::ComputeApproximateAngle_3 */

