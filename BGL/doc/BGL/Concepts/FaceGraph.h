/*!
\ingroup PkgBGLConcepts
\cgalConcept

The concept `FaceGraph` refines the concept `HalfedgeGraph`. 
It adds the requirements for a graph to explicitly
maintain faces described by halfedges, to provide access from a face to
an incident halfedge, and to provide access from a halfedge to its incident
face. 

\cgalRefines `HalfedgeGraph`
\cgalHasModel `CGAL::Polyhedron_3`
\cgalHasModel `CGAL::Surface_mesh`
\cgalHasModel `CGAL::Linear_cell_complex_for_combinatorial_map`

*/
class FaceGraph {};

/*! \relates FaceGraph
returns the face incident to halfedge `h`.
 */
template <typename FaceGraph>
boost::graph_traits<FaceGraph>::face_descriptor
face(boost::graph_traits<FaceGraph>::halfedge_descriptor h, const FaceGraph& g);

/*! \relates FaceGraph
returns the halfedge incident to face `f`.
 */
template <typename FaceGraph>
boost::graph_traits<FaceGraph>::halfedge_descriptor
halfedge(boost::graph_traits<FaceGraph>::face_descriptor f, const FaceGraph& g);

/*! \relates FaceGraph
returns the number of halfedges incident to face `f`.
 */
template <typename FaceGraph>
boost::graph_traits<FaceGraph>::degree_size_type
degree(boost::graph_traits<FaceGraph>::face_descriptor f, const FaceGraph& g);

/*! \relates FaceGraph
returns a special face that is not equal to any other face.
 */
template <typename FaceGraph>
boost::graph_traits<FaceGraph>::face_descriptor
null_face(const FaceGraph& g);
