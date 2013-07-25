
/*!
\ingroup PkgSphericalKernel3GeometricConcepts
\cgalConcept

\brief

\cgalRefines `Kernel::HasOnUnboundedSide_3`

\sa `SphericalKernel::BoundedSide_3`
\sa `SphericalKernel::HasOnBoundedSide_3`

*/

class SphericalKernel::HasOnUnboundedSide_3 {
public:

/// \name Operations 
/// An object of this type must provide:
/// @{

/*!
For a sphere. 
*/ 
bool 
operator() 
(const SphericalKernel::Sphere_3& s, 
const SphericalKernel::Circular_arc_point_3& p); 

/*!
For a circle. 
*/ 
bool 
operator() 
(const SphericalKernel::Circle_3& s, 
const SphericalKernel::Circular_arc_point_3& p); 

/// @}

}; /* end SphericalKernel::HasOnUnboundedSide_3 */

