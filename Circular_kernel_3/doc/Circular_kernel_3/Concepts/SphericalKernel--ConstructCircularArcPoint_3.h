
/*!
\ingroup PkgSphericalKernel3GeometricConcepts
\cgalConcept

*/

class SphericalKernel::ConstructCircularArcPoint_3 {
public:

/// \name Operations
/// A model of this concept must provide: 
/// @{

/*!

*/ 
SphericalKernel::Circular_arc_point_3 operator() 
(const SphericalKernel::Root_for_spheres_2_3 & r); 

/*!

*/ 
SphericalKernel::Circular_arc_point_3 operator() 
(const SphericalKernel::Point_3 & p); 

/// @}

}; /* end SphericalKernel::ConstructCircularArcPoint_3 */

