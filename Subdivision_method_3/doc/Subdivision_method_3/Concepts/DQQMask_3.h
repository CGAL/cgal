
/*!
\ingroup PkgSurfaceSubdivisionMethod3Concepts
\cgalConcept

Required member functions for the `DQQMask_3` concept. This
policy concept of geometric computations is used in
`CGAL::Subdivision_method_3::DQQ<PolygonMesh, Mask, NamedParameters>`.

\image html DSCornerMask.svg

\cgalRefines `SubdivisionMask_3`

\cgalHasModel `CGAL::DooSabin_mask_3<PolygonMesh, VertexPointMap>`

\sa `CGAL::Subdivision_method_3`

*/
class DQQMask_3 {
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
DQQMask_3(PolygonMesh* pmesh);

/*! Constructor.
*/
DQQMask_3(PolygonMesh* pmesh, VertexPointMap vpmap);

/*!
computes the subdivided point `pt` based on the neighborhood
of the vertex pointed by the halfedge `hd`.
*/
void corner_node(halfedge_descriptor hd, Point& pt);

/// @}

}; /* end DQQMask_3 */

