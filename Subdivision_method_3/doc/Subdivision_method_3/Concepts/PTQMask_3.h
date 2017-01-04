/*!
\ingroup PkgSurfaceSubdivisionMethods3Concepts
\cgalConcept

Required member functions for the `PTQMask_3` concept. This 
policy concept of geometric computations is used in 
`CGAL::Subdivision_method_3::PTQ<PolygonMesh, Mask>`. 

\image html LoopBorderMask.png
\image latex LoopBorderMask.png

\cgalRefines `SubdivisionMask_3`

\cgalHasModel `CGAL::Loop_mask_3<PolygonMesh>`

\sa `CGAL::Subdivision_method_3`
*/
class PTQMask_3 {
public:

/// \name Types 
/// @{

/*!  The polygon mesh must be triangulated.

*/ 
  typedef unspecified_type PolygonMesh;

/// @}


/// \name Operations 
/// @{

/*! Constructor
 */

PTQMask_3(PolygonMesh& pmesh);

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

}; /* end PTQMask_3 */

