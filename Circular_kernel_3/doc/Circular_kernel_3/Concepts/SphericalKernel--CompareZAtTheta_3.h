
/*!
\ingroup PkgSphericalKernel3Concepts
\cgalconcept

\sa `SphericalKernel::CompareZToRight_3`

*/
class SphericalKernel::CompareZAtTheta_3 {
public:

/// \name Operations
/// An object `fo` of this type must provide: 
/// @{

/*! 

compares the \f$ z\f$-coordinates of the two intersections points of \f$ a0\f$ and \f$ a1\f$ with the meridian defined by \f$ m\f$ (see section \ref sectionSKobjects). 
\pre `a0` and `a1` lie on the context sphere used by the function `SphericalKernel::compare_z_at_theta_3_object`. \f$ m \neq(0,0,0)\f$ and the \f$ z\f$-coordinate of \f$ m\f$ is \f$ 0\f$. Arcs \f$ a0\f$ and \f$ a1\f$ are \f$ \theta\f$-monotone and both intersected by the meridian defined by \f$ m\f$ (see section \ref sectionSKobjects). 
*/ 
Comparison_result operator() 
( const SphericalKernel::Circular_arc_3& a0, 
const SphericalKernel::Circular_arc_3& a1, 
const SphericalKernel::Vector_3 &m); 

/*! 
given a meridian anchored at the poles of the context sphere used by the function `SphericalKernel::compare_z_at_theta_3_object`, and passing through point \f$ p\f$, 
compares the \f$ z\f$-coordinate of point \f$ p\f$ and that of the intersection of the meridian with \f$ a\f$. 
\pre `a` and `p` lie on the context sphere used by the function `SphericalKernel::compare_z_at_theta_3_object`, arc \f$ a\f$ is \f$ \theta\f$-monotone and the meridian passing through \f$ p\f$ intersects arc \f$ a\f$. 
*/ 
Comparison_result operator() 
( const SphericalKernel::Circular_arc_point_3& p, 
const SphericalKernel::Circular_arc_3& a); 

/// @}

}; /* end SphericalKernel::CompareZAtTheta_3 */

