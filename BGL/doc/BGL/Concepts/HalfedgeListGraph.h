/*!
\ingroup PkgBGLConcepts
\cgalConcept

The concept `HalfedgeListGraph` refines the concept `HalfedgeGraph`
and adds the requirements for traversal of all halfedges in the graph.

\cgalAssociatedTypesBegin

\cgalAssociatedTypeBegin{boost::graph_traits<HalfedgeListGraph>::%halfedge_iterator} A halfedge iterator (obtained via `halfedges(g)`) provides access to all of the halfedges in a graph.
A halfedge iterator type must meet the requirements of `MultiPassInputIterator`. The value type of the
halfedge iterator must be the same as the halfedge descriptor of the graph.
\cgalAssociatedTypeEnd

\cgalAssociatedTypesEnd

\cgalRefines `HalfedgeGraph`

\cgalHasModel See \link PkgBGLTraits Boost Graph Traits Specializations \endlink

\sa \link PkgBGLConcepts Graph Concepts \endlink
*/

class HalfedgeListGraph {};

/*! \relates HalfedgeListGraph
 * returns an iterator range over all halfedges.
 */
template <typename HalfedgeListGraph>
std::pair<boost::graph_traits<HalfedgeListGraph>::halfedge_iterator,
          boost::graph_traits<HalfedgeListGraph>::halfedge_iterator>
halfedges(const HalfedgeListGraph& g);


/*! \relates HalfedgeListGraph
  returns an upper bound of the number of halfedges of the graph.
  \attention `num_halfedges()` may return a number larger than `std::distance(halfedges(g).first, halfedges(g).second)`.
  This is the case for implementations only marking halfedges deleted in the halfedge container.
 */
template <typename HalfedgeListGraph>
boost::graph_traits<HalfedgeListGraph>::halfedge_size_type
num_halfedges(const HalfedgeListGraph& g);

