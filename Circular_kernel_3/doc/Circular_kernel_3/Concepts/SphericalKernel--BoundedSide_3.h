
/*!
\ingroup PkgSphericalKernel3GeometricConcepts
\cgalConcept

\brief

\cgalRefines `Kernel::BoundedSide_3`

\sa `SphericalKernel::HasOnBoundedSide_3`
\sa `SphericalKernel::HasOnUnboundedSide_3`

*/

class SphericalKernel::BoundedSide_3 {
public:

/// \name Operations 
/// An object of this type must provide:
/// @{

/*!
For a sphere. 
*/ 
Bounded_side 
operator() 
(const SphericalKernel::Sphere_3& s, 
const SphericalKernel::Circular_arc_point_3& p); 

/*!
For a circle. 
*/ 
Bounded_side 
operator() 
(const SphericalKernel::Circle_3& s, 
const SphericalKernel::Circular_arc_point_3& p); 

/// @}

}; /* end SphericalKernel::BoundedSide_3 */

