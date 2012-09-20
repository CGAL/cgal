
/*!
\ingroup PkgSphericalKernel3GeometricConcepts
\cgalconcept

\brief Testing whether two curves intersect. 

\refines `Kernel::DoIntersect_3`

\sa `Kernel::do_intersect`"CGAL::do_intersect" 

*/

class SphericalKernel::DoIntersect_3 {
public:

/// \name Operations 
/// An object `fo` of this type must provide: 
/// @{

/*! 
determines if two geometric objects of type Type1 and Type2 intersect or not. 

for all pairs `Type1` and `Type2`, where the types `Type1` and `Type2`
can be either, any of the following: 

- `SphericalKernel::Plane_3`
- `SphericalKernel::Line_3`
- `SphericalKernel::Line_arc_3`
- `SphericalKernel::Sphere_3`
- `SphericalKernel::Circle_3`
*/ 
bool operator() 
(const Type1 & obj1, const Type2 & obj2); 

/// @}

}; /* end SphericalKernel::DoIntersect_3 */

