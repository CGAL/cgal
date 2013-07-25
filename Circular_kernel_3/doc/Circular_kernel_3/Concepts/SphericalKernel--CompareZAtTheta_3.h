
/*!
\ingroup PkgSphericalKernel3GeometricConcepts
\cgalConcept

\sa `SphericalKernel::CompareZToRight_3`

*/
class SphericalKernel::CompareZAtTheta_3 {
public:

/// \name Operations
/// An object of this type must provide: 
/// @{

/*!

compares the \f$ z\f$-coordinates of the two intersections points of `a0` and `a1` with the meridian defined by `m` (see Section \ref sectionSKobjects). 
\pre `a0` and `a1` lie on the context sphere used by the function `SphericalKernel::compare_z_at_theta_3_object`. `m` \f$ \neq(0,0,0)\f$ and the \f$ z\f$-coordinate of `m` is \f$ 0\f$. Arcs `a0` and `a1` are \f$ \theta\f$-monotone and both intersected by the meridian defined by `m`(see Section \ref sectionSKobjects). 
*/ 
Comparison_result operator() 
( const SphericalKernel::Circular_arc_3& a0, 
const SphericalKernel::Circular_arc_3& a1, 
const SphericalKernel::Vector_3 &m); 

/*!
given a meridian anchored at the poles of the context sphere used by the function `SphericalKernel::compare_z_at_theta_3_object`, and passing through point `p`, 
compares the \f$ z\f$-coordinate of point `p` and that of the intersection of the meridian with `a`. 
\pre `a` and `p` lie on the context sphere used by the function `SphericalKernel::compare_z_at_theta_3_object`, arc `a` is \f$ \theta\f$-monotone and the meridian passing through `p` intersects arc `a`. 
*/ 
Comparison_result operator() 
( const SphericalKernel::Circular_arc_point_3& p, 
const SphericalKernel::Circular_arc_3& a); 

/// @}

}; /* end SphericalKernel::CompareZAtTheta_3 */

