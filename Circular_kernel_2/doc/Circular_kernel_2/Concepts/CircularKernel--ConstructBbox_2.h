
/*!
\ingroup PkgCircularKernel2GeometricConcepts
\cgalconcept

A model `fo` of this type must provide operators to construct 
a bounding box of geometric objects: 

*/

class CircularKernel::ConstructBbox_2 {
public:

/// \name See Also 
/// @{

/*! 

*/ 
CGAL::Bbox_2 operator() 
(const CircularKernel::Circular_arc_point_2 & p); 

/*! 

*/ 
CGAL::Bbox_2 operator() 
(const CircularKernel::Line_arc_2 & l); 

/*! 

*/ 
CGAL::Bbox_2 operator() 
(const CircularKernel::Circular_arc_2 & c); 

/// @}

}; /* end CircularKernel::ConstructBbox_2 */

