
/*!
\ingroup PkgSphericalKernel3GeometricConcepts
\cgalConcept

\cgalRefines Kernel::ConstructBbox_3

*/

class SphericalKernel::ConstructBbox_3 {
public:

/// \name Operations
/// A model of this concept must provide operators to construct 
/// a bounding box of geometric objects: 
/// @{

/*!

*/ 
CGAL::Bbox_3 operator() 
(const SphericalKernel::Circular_arc_point_3 & p); 

/*!

*/ 
CGAL::Bbox_3 operator() 
(const SphericalKernel::Line_arc_3 & l); 

/*!

*/ 
CGAL::Bbox_3 operator() 
(const SphericalKernel::Circular_arc_3 & a); 

/// @}

}; /* end SphericalKernel::ConstructBbox_3 */

