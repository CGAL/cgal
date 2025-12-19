/*!
\ingroup PkgBGLConcepts
\cgalConcept

Concept from the Boost Graph Library.
See https://www.boost.org/libs/graph/doc/VertexListGraph.html.

The concept `VertexListGraph` refines the concept
 <a href="https://www.boost.org/libs/graph/doc/Graph.html"><code>Graph</code></a>
and adds the requirement for traversal of all vertices in a graph.

\cgalRefines{<a href="https://www.boost.org/libs/graph/doc/Graph.html">Graph</a>}

\cgalHasModelsBegin
\cgalHasModelsBare{See \link PkgBGLTraits Boost Graph Traits Specializations \endlink}
\cgalHasModelsEnd

\sa \link PkgBGLConcepts Graph Concepts \endlink
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
boost::graph_traits<VertexListGraph>::vertices_size_type
num_vertices(const VertexListGraph& g);

