
/*!
\ingroup PkgSurfaceSubdivisionMethod3Concepts
\cgalConcept

Required member functions for the `Sqrt3Mask_3` concept. This
policy concept of geometric computations is used in
`CGAL::Subdivision_method_3::Sqrt3<PolygonMesh, Mask, NamedParameters>`.

\cgalRefines `SubdivisionMask_3`

\cgalHasModel `CGAL::Sqrt3_mask_3<PolygonMesh, VertexPointMap>`

\sa `CGAL::Subdivision_method_3`

*/

class Sqrt3Mask_3 {
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

/*! Constructor.
 * The default vertex point property map is used.
 */
Sqrt3Mask_3(PolygonMesh* pmesh);

/*! Constructor
 */
Sqrt3Mask_3(PolygonMesh* pmesh, VertexPointMap vpmap);

/*!
computes the subdivided point `pt` based on the neighborhood
of the face `fd`.
*/
void face_node(face_descriptor fd, Point_3& pt);

/*!
computes the subdivided point `pt` based on the neighborhood
of the vertex `vd`.
*/
void vertex_node(vertex_descriptor vd, Point& pt);

/*!
computes the subdivided points `ept1` and `ept2` based on the neighborhood
of the halfedge `hd` (whose opposite is on the border).
Along `hd`, `ept1` comes before `ept2`.
`vpt` is the updated point for the target vertex of `hd`.
*/
void border_node(halfedge_descriptor hd, Point& ept1, Point& ept2, Point& vpt);

/// @}

}; /* end Sqrt3Mask_3 */

