/*!
\ingroup PkgSurfaceSubdivisionMethods3Concepts
\cgalConcept

Required member functions for the `PTQMask_3` concept. This 
policy concept of geometric computations is used in 
`CGAL::Subdivision_method_3::PTQ<Polyhedron_3, Mask>`. 

\image html LoopBorderMask.png
\image latex LoopBorderMask.png

\cgalHasModel `CGAL::Loop_mask_3<Polyhedron_3>`

\sa `CGAL::Subdivision_method_3`
*/
class PTQMask_3 {
public:

/// \name Operations 
/// @{

/*!
computes the edge-point `pt` based on the neighborhood 
of the edge `e`. 
*/ 
void edge_node(Edge_handle e, Point_3& pt); 

/*!

computes the vertex-point `pt` based on the neighborhood 
of the vertex `v`. 
*/ 
void vertex_node(Vertex_handle v, Point_3& pt); 

/*!

computes the edge-point `ept` and the vertex-point `vpt` 
based on the neighborhood of the border edge `e`. 
*/ 
void border_node(Halfedge_handle e, Point_3& ept, Point_3& vpt); 

/// @}

}; /* end PTQMask_3 */

