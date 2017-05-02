/*!
\ingroup PkgBGLConcepts
\cgalConcept

The concept `EdgeListGraph` refines the concept `Graph` and adds
the requirement for traversal of all edges in a graph.

\cgalRefines `Graph`
\cgalHasModel `CGAL::Polyhedron_3`
\cgalHasModel `CGAL::Surface_mesh`

*/
class EdgeListGraph{};


/*! \relates EdgeListGraph
 * returns an iterator range over all edges.
 */
std::pair<boost::graph_traits<EdgeListGraph>::edge_iterator,
          boost::graph_traits<EdgeListGraph>::edge_iterator>
edges(const EdgeListGraph& g);


/*! \relates EdgeListGraph
  returns an upper bound of the number of edges of the graph.
  \attention `num_edges()` may return a number larger than `std::distance(edges(g).first, edges(g).second)`.
  This is the case for implementations only marking edges deleted in the edge container.
 */
boost::graph_traits<EdgeListGraph>::ver_size_type
num_edges(const EdgeListGraph& g);


/*! \relates EdgeGraph
returns the source vertex of `h`.
 */
boost::graph_traits<EdgeGraph>::vertex_descriptor
source(boost::graph_traits<EdgeGraph>::halfedge_descriptor h, EdgeGraph& g);


/*! \relates EdgeGraph
returns the target vertex of `h`.
 */
boost::graph_traits<EdgeGraph>::vertex_descriptor
target(boost::graph_traits<EdgeGraph>::halfedge_descriptor h, EdgeGraph& g);
