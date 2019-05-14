
/*!
\ingroup PkgSphericalKernel3GeometricConcepts
\cgalConcept


\sa `SphericalKernel::CompareZAtTheta_3`

*/

class SphericalKernel::CompareZToRight_3 {
public:

/// \name Operations
/// An object of this type must provide: 
/// @{

/*!
Compares the \f$ z\f$-coordinates of the intersection points of both arcs 
with a meridian anchored at the poles of the context sphere used by the function `SphericalKernel::compare_z_to_right_3_object`, at a \f$ \theta\f$-coordinate 
infinitesimally greater that the \f$ \theta\f$-coordinate of point `p`. 
\pre `a0` and `a1` lie on the context sphere used by the function `SphericalKernel::compare_z_to_right_3_object`, `a0` and `a1` are \f$ \theta\f$-monotone, `p` lies on `a0` and `a1` and is not a \f$ \theta\f$-extremal point of the supporting circle of `a0` or `a1`. 
*/ 
Comparison_result operator() 
( const SphericalKernel::Circular_arc_3& a0, 
const SphericalKernel::Circular_arc_3& a1, 
const SphericalKernel::Circular_arc_point_3 &p); 

/// @}

}; /* end SphericalKernel::CompareZToRight_3 */

