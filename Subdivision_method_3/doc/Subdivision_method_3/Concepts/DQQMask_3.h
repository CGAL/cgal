
/*!
\ingroup PkgSurfaceSubdivisionMethods3Concepts
\cgalConcept

Required member functions for the `DQQMask_3` concept. This 
policy concept of geometric computations is used in 
`CGAL::Subdivision_method_3::DQQ<Polyhedron_3, Mask>`. 

\image html DSCornerMask.png
\image latex DSCornerMask.png

\cgalRefines `SubdivisionMask`

\cgalHasModel `CGAL::DooSabin_mask_3<Polyhedron_3>`

\sa `CGAL::Subdivision_method_3`

*/
class DQQMask_3 {
public:

/// \name Operations 
/// @{

/*! Constructor
 */

  DQQMask_3(PolygonMesh& pm);

/*!
computes the subdivided point `pt` based on the neighborhood 
of the vertex pointed by the halfedge `he`. 
*/ 
void corner_node(halfedge_descriptor he, Point& pt); 

/// @}

}; /* end DQQMask_3 */

