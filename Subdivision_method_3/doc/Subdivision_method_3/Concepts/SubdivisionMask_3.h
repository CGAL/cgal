
/*!
\ingroup PkgSurfaceSubdivisionMethod3Concepts
\cgalConcept

This concept defines types used by the mask concepts.

\sa `CGAL::Subdivision_method_3`

*/

class SubdivisionMask_3 {
public:

/// \name Types
/// @{

/*!

*/
  typedef unspecified_type PolygonMesh;

  typedef boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
  typedef boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  typedef boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;

/// @}

}; /* end PQQMask_3 */

