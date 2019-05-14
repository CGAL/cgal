
/*!
\ingroup PkgSphericalKernel3GeometricConcepts
\cgalConcept

*/
class SphericalKernel::ConstructCircularSourceVertex_3 {
public:

/// \name Operations
/// A model of this concept must provide: 
/// @{

/*!
Constructs the source vertex of `a`. 
*/ 
SphericalKernel::Circular_arc_point_3 operator() 
(const SphericalKernel::Circular_arc_3 & a); 

/*!
Same, for a line segment. 
*/ 
SphericalKernel::Circular_arc_point_3 operator() 
(const SphericalKernel::Line_arc_3 & l); 

/// @}

}; /* end SphericalKernel::ConstructCircularSourceVertex_3 */

