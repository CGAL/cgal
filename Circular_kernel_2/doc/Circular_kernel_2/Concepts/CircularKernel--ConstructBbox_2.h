
/*!
\ingroup PkgCircularKernel2GeometricConcepts
\cgalConcept

A function object concept to construct a bounding box of geometric objects:

\cgalRefines Kernel::ConstructBbox_2

*/

class CircularKernel::ConstructBbox_2 {
public:

/// \name Operations
/// A model of this concept must provide:
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

