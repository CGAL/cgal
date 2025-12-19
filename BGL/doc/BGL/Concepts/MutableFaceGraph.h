/*!
\ingroup PkgBGLConcepts
\cgalConcept

The concept `MutableFaceGraph` refines the concepts `FaceGraph` and `MutableHalfedgeGraph` and adds
the requirement for operations to add faces and to modify face-halfedge relations.

\cgalRefines{FaceGraph,MutableHalfedgeGraph}

\cgalHasModelsBegin
\cgalHasModelsBare{See \link PkgBGLTraits Boost Graph Traits Specializations \endlink}
\cgalHasModelsEnd

\sa \link PkgBGLConcepts Graph Concepts \endlink
*/
class MutableFaceGraph{};

/*! \relates MutableFaceGraph
Adds a new face to the graph without initializing the connectivity.
 */
template <typename MutableFaceGraph>
boost::graph_traits<MutableFaceGraph>::face_descriptor
add_face(MutableFaceGraph& g);

/*! \relates MutableFaceGraph
Removes `f` from the graph.
 */
template <typename MutableFaceGraph>
boost::graph_traits<MutableFaceGraph>::face_descriptor
remove_face(boost::graph_traits<MutableFaceGraph>::face_descriptor f, MutableFaceGraph& g);

/*! \relates MutableFaceGraph
Sets the corresponding face of `h` to `f`.
 */
template <typename MutableFaceGraph>
void
set_face(boost::graph_traits<MutableFaceGraph>::halfedge_descriptor h, boost::graph_traits<MutableFaceGraph>::face_descriptor f, MutableFaceGraph& g);

/*! \relates MutableFaceGraph
Sets the corresponding halfedge of `f` to `h`.
 */
template <typename MutableFaceGraph>
void
set_halfedge(boost::graph_traits<MutableFaceGraph>::face_descriptor f, boost::graph_traits<MutableFaceGraph>::halfedge_descriptor h, MutableFaceGraph& g);

/*! \relates MutableFaceGraph
Indicates the expected size of vertices (`nv`), edges (`ed`) and faces (`nf`).
 */
template <typename MutableFaceGraph>
void
reserve(MutableFaceGraph& g, boost::graph_traits<MutableFaceGraph>::vertices_size_type nv, boost::graph_traits<MutableFaceGraph>::edges_size_type ne, boost::graph_traits<MutableFaceGraph>::faces_size_type nf);
