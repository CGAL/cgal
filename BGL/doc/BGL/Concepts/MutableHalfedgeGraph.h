/*!
\ingroup PkgBGLConcepts
\cgalConcept

The concept `MutableHalfedgeGraph` refines the concept `HalfedgeGraph`
and adds the requirements for operations to add vertices and edges, and to
update the incidence information between vertices and halfedges.

\cgalRefines `HalfedgeGraph`

\cgalHasModel See \link PkgBGLTraits Boost Graph Traits Specializations \endlink

\sa \link PkgBGLConcepts Graph Concepts \endlink
*/
class MutableHalfedgeGraph{};


/*! \relates MutableFaceGraph
Adds a new vertex to the graph without initializing the connectivity.
 */
template <typename MutableHalfedgeGraph>
boost::graph_traits<MutableHalfedgeGraph>::face_descriptor
add_vertex(MutableHalfedgeGraph& g);

/*! \relates MutableHalfedgeGraph
Removes `v` from the graph.
 */
template <typename MutableHalfedgeGraph>
boost::graph_traits<MutableHalfedgeGraph>::face_descriptor
remove_vertex(boost::graph_traits<MutableHalfedgeGraph>::vertex_descriptor v, MutableHalfedgeGraph& g);

/*! \relates MutableFaceGraph
Adds two opposite halfedges to the graph without initializing the connectivity.
 */
template <typename MutableHalfedgeGraph>
boost::graph_traits<MutableHalfedgeGraph>::edge_descriptor
add_edge(MutableHalfedgeGraph& g);

/*! \relates MutableHalfedgeGraph
Removes the two halfedges corresponding to `e` from the graph.
 */
template <typename MutableHalfedgeGraph>
boost::graph_traits<MutableHalfedgeGraph>::face_descriptor
remove_edge(boost::graph_traits<MutableHalfedgeGraph>::edge_descriptor e, MutableHalfedgeGraph& g);


/*! \relates MutableHalfedgeGraph
Sets the target vertex of `h` and the source of `opposite(h)` to `v`.
 */
template <typename MutableHalfedgeGraph>
void
set_target(boost::graph_traits<MutableHalfedgeGraph>::halfedge_descriptor h, boost::graph_traits<MutableHalfedgeGraph>::vertex_descriptor v, MutableHalfedgeGraph& g);

/*! \relates MutableHalfedgeGraph
Sets the halfedge of `v` to `h`. The target vertex of `h` must be `v`.
 */
template <typename MutableHalfedgeGraph>
void
set_halfedge(boost::graph_traits<MutableHalfedgeGraph>::vertex_descriptor v, boost::graph_traits<MutableHalfedgeGraph>::halfedge_descriptor h, MutableHalfedgeGraph& g);


/*! \relates MutableHalfedgeGraph
Sets the successor of `h1` around a face to `h2`, and the prededecessor of `h2` to `h1`.
 */
template <typename MutableHalfedgeGraph>
void
set_next(boost::graph_traits<MutableHalfedgeGraph>::halfedge_descriptor h1, boost::graph_traits<MutableHalfedgeGraph>::halfede_descriptor h2, MutableHalfedgeGraph& g);
