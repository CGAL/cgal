/*!
\ingroup PkgBGLConcepts
\cgalConcept

The concept `VertexListGraph` refines the concept `Graph` and adds
the requirement for traversal of all vertices in a graph.

\cgalRefines `Graph`
\cgalHasModel `CGAL::Polyhedron_3`
\cgalHasModel `CGAL::Surface_mesh`

*/
class VertexListGraph{};


/*! \relates VertexListGraph
 * returns an iterator range over all vertices.
 */
template <typename VertexListGraph>
std::pair<boost::graph_traits<VertexListGraph>::vertex_iterator,
          boost::graph_traits<VertexListGraph>::vertex_iterator>
vertices(const VertexListGraph& g);


/*! \relates VertexListGraph
  returns an upper bound of the number of vertices of the graph.
  \attention `num_vertices()` may return a number larger than `std::distance(vertices(g).first, vertices(g).second)`.
  This is the case for implementations only marking vertices deleted in the vertex container.
 */
template <typename VertexListGraph>
boost::graph_traits<VertexListGraph>::ver_size_type
num_vertices(const VertexListGraph& g);

