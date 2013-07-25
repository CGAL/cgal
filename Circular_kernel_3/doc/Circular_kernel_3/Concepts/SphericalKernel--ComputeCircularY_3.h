
/*!
\ingroup PkgSphericalKernel3GeometricConcepts
\cgalConcept
*/
class SphericalKernel::ComputeCircularY_3 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
returns the \f$ y\f$-coordinate of the point `p`. 
*/ 

SphericalKernel::Root_of_2 
operator()(const SphericalKernel::Circular_arc_point_3 & p) const; 

///@}

}; /* end SphericalKernel::ComputeCircularY_3 */

