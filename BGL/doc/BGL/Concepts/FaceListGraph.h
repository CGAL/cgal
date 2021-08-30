/*!
\ingroup PkgBGLConcepts
\cgalConcept

The concept `FaceListGraph` refines the concept `FaceGraph` and adds
the requirement for traversal of all faces in a graph.

\cgalAssociatedTypesBegin

\cgalAssociatedTypeBegin{boost::graph_traits<FaceListGraph>::%face_iterator} A face iterator (obtained via `faces(g)`) provides access to all of the faces in a graph.
A face iterator type must meet the requirements of `MultiPassInputIterator`. The value type of the
face iterator must be the same as the face descriptor of the graph.
\cgalAssociatedTypeEnd

\cgalAssociatedTypesEnd

\cgalRefines `FaceGraph`

\cgalHasModel See \link PkgBGLTraits Boost Graph Traits Specializations \endlink

\sa \link PkgBGLConcepts Graph Concepts \endlink
*/
class FaceListGraph{};


/*! \relates FaceListGraph
 * returns an iterator range over all faces.
 */
template <typename FaceListGraph>
std::pair<boost::graph_traits<FaceListGraph>::face_iterator,
          boost::graph_traits<FaceListGraph>::face_iterator>
faces(const FaceListGraph& g);


/*! \relates FaceListGraph
  returns an upper bound of the number of faces of the graph.
  \attention `num_faces()` may return a number larger than `std::distance(faces(g).first, faces(g).second)`.
  This is the case for implementations only marking faces deleted in the face container.
 */
template <typename FaceListGraph>
boost::graph_traits<FaceListGraph>::face_size_type
num_faces(const FaceListGraph& g);

