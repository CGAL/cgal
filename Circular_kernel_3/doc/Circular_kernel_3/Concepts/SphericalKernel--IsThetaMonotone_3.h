
/*!
\ingroup PkgSphericalKernel3GeometricConcepts
\cgalConcept


\sa `SphericalKernel::MakeThetaMonotone_3`

*/

class SphericalKernel::IsThetaMonotone_3 {
public:

/// \name Operations
/// An object of this type must provide: 
/// @{

/*!
Tests whether the arc `a` is \f$ \theta\f$-monotone, i.e.\ the intersection of 
any meridian anchored at the poles of the context sphere used by the function `SphericalKernel::is_theta_monotone_3_object` 
and the arc `a` is reduced to at most one point in general, and two points if a pole of that sphere is 
an endpoint of `a`. Note that a bipolar circle has no such arcs. 
\pre `a` lies on the context sphere used by the function `SphericalKernel::is_theta_monotone_3_object`, and the supporting circle of `a` is not bipolar. 

*/ 
bool operator() 
(const SphericalKernel::Circular_arc_3 &a); 

/// @}

}; /* end SphericalKernel::IsThetaMonotone_3 */

