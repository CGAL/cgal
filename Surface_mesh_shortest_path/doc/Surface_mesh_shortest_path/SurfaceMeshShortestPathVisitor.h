/*!
\ingroup PkgSurfaceMeshShortestPathConcepts

\cgalConcept

The concept `SurfaceMeshShortestPathVisitor` describes the visitor type
used to collect the edges and vertices traversed by a shortest path on the
surface of a polyhedron. 

The order of the method calls reports the sequence of the shortest path from 
the query point back to the nearest source point, excluding the query point itself.

*/

class SurfaceMeshShortestPathVisitor
{
public:

/// \name Methods
/// @{

  /*!
  \brief Called when an edge was traversed in the shortest path sequence.
  \param edge halfedge of the polyhedron crossed by the shortest path. The halfedge given will be directed to the face <em>nearest</em> the query point.
  \param t value in the range [0.0,1.0] specifying where along `edge` the shortest path crossed. 
    - 0.0 means it crossed at `source(edge)`
    - 1.0 means it crossed at `target(edge)`
    - Any other value is linearly interpolated between the endpoints.
    Note that values of 0.0 and 1.0 are possible in some degenerate situations.
  */
  void edge(halfedge_descriptor edge, FT t);
  
  /*!
  \brief Called when a vertex is encountered in the shortest path sequence.
  \param vertex vertex of the polyhedron encountered by the shortest path
  */
  void vertex(vertex_descriptor vertex);
  
  /*!
  \brief Called when a face location is encountered in the shortest path sequence.
  \remarks This will only be called as the last element in the path sequence, and only
    if the source point is an internal face location (i.e. not an edge or a vertex).
  \param face face of the polyhedron encountered by the shortest path
  \param location the barycentric coordinate inside face of this point on the path
  */
  void face(face_descriptor face, Barycentric_coordinate location);
};