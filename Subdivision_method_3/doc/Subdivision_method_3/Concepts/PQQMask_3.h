
/*!
\ingroup PkgSurfaceSubdivisionMethods3Concepts
\cgalConcept

Required member functions for the `PQQMask_3` concept. This 
policy concept of geometric computations is used in 
`CGAL::Subdivision_method_3::PQQ<PolygonMesh, Mask>`. 

\image html CCBorderMask.png
\image latex CCBorderMask.png

\cgalRefines `SubdivisionMask_3`

\cgalHasModel `CGAL::CatmullClark_mask_3<PolygonMesh>`

\sa `CGAL::Subdivision_method_3`

*/

class PQQMask_3 {
public:

/// \name Operations 
/// @{

/*! Constructor
 */

  PQQMask_3(PolygonMesh& pmesh);

/*!
computes the face-point `pt` based on the neighborhood 
of the face `fd`. 
*/ 
void face_node(face_descriptor fd, Point_3& pt); 

/*!

computes the edge-point `pt` based on the neighborhood 
of the edge `hd`. 
*/ 
void edge_node(halfedge_descriptor hd, Point_3& pt); 

/*!

computes the vertex-point `pt` based on the neighborhood 
of the vertex `vd`. 
*/ 
void vertex_node(vertex_descriptor vd, Point_3& pt); 

/*!

computes the edge-point `ept` and the vertex-point `vpt` 
based on the neighborhood of the border edge `hd`. 
*/ 
void border_node(halfedge_descriptor hd, Point_3& ept, Point_3& vpt); 

/// @}

}; /* end PQQMask_3 */

