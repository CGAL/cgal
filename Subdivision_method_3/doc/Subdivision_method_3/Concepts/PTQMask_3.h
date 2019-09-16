/*!
\ingroup PkgSurfaceSubdivisionMethod3Concepts
\cgalConcept

Required member functions for the `PTQMask_3` concept. This
policy concept of geometric computations is used in
`CGAL::Subdivision_method_3::PTQ<PolygonMesh, Mask, NamedParameters>`.

\image html LoopBorderMask.svg

\cgalRefines `SubdivisionMask_3`

\cgalHasModel `CGAL::Loop_mask_3<PolygonMesh, VertexPointMap>`

\sa `CGAL::Subdivision_method_3`
*/
class PTQMask_3 {
public:

/// \name Types
/// @{

/*!  The polygon mesh must be triangulated.

*/
  typedef unspecified_type PolygonMesh;
  typedef unspecified_type VertexPointMap;
/// @}


/// \name Operations
/// @{

/*! Constructor
 *
 * The default vertex point property map is used.
*/
PTQMask_3(PolygonMesh* pmesh);

/*! Constructor
*/
PTQMask_3(PolygonMesh* pmesh, VertexPointMap vpmap);

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
based on the neighborhood of the border edge of `hd`. `hd` is not
a border halfedge (its opposite is) and `vpt` corresponds to the
target vertex of `hd`.
*/
void border_node(halfedge_descriptor hd, Point_3& ept, Point_3& vpt);

/// @}

}; /* end PTQMask_3 */

