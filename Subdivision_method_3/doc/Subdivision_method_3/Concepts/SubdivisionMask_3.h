
/*!
\ingroup PkgSurfaceSubdivisionMethods3Concepts
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
  typedef unspecified_type Mesh;
  typedef Mesh TriangleMesh;
  typedef Mesh PolygonMesh;
  typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
  typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
  typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
  

/// @}

}; /* end PQQMask_3 */

