/*!
\ingroup PkgBGLConcepts
\cgalConcept

The concept `HalfedgeGraph` is a refinement of the \bgl concept
`IncidenceGraph` and adds the notion of a *halfedge*: Each edge is
associated with two *opposite* halfedges with source and target vertices swapped.
Furthermore, halfedges have a *successor* and *predecessor*,
and form cycles we call *faces*. However, this concept
does not introduce a face type.
A `HalfedgeGraph` is undirected and does not allow parallel edges.

Using the composition of the *successor* and *opposite* functions results
in another cycle, namely the cycle of halfedges which are incident to
the same vertex. We refer to \ref PkgBGLIterators for a description of
iterators and circulators for these halfedge cycles.

A partial specialization must be provided for `CGAL::graph_has_property`
for each internal property map available.

\cgalAssociatedTypesBegin

\cgalAssociatedTypeBegin{boost::graph_traits<HalfedgeGraph>::%vertex_descriptor} A vertex descriptor corresponds to a unique vertex in an abstract graph instance.
A vertex descriptor must be `DefaultConstructible`, `Assignable`, `EqualityComparable`, and `Hashable`.
\cgalAssociatedTypeEnd

\cgalAssociatedTypeBegin{boost::graph_traits<HalfedgeGraph>::%halfedge_descriptor} A halfedge descriptor corresponds to a unique halfedge in an abstract graph instance.
A halfedge descriptor must be `DefaultConstructible`, `Assignable`, `EqualityComparable`, and `Hashable`.
\cgalAssociatedTypeEnd

\cgalAssociatedTypeBegin{boost::graph_traits<HalfedgeGraph>::%edge_descriptor} An edge descriptor corresponds to a unique edge in an abstract graph instance.
An edge descriptor must be `DefaultConstructible`, `Assignable`, `EqualityComparable`, and `Hashable`.
\cgalAssociatedTypeEnd

\cgalAssociatedTypesEnd

\cgalRefines{IncidenceGraph,PropertyGraph}

A model of `HalfedgeGraph` must have the interior property `vertex_point` attached to its vertices.

\cgalHasModelsBegin
\cgalHasModelsBare{See \link PkgBGLTraits Boost Graph Traits Specializations \endlink}
\cgalHasModelsEnd

\sa \link PkgBGLConcepts Graph Concepts \endlink
*/
class HalfedgeGraph {
  /// Returns a special `boost::graph_traits<HalfedgeGraph>::%halfedge_descriptor` object which
  /// does not refer to any halfedge of graph object which type is `HalfedgeGraph`.
  static boost::graph_traits<HalfedgeGraph>::halfedge_descriptor null_halfedge();
};

/*! \relates HalfedgeGraph
returns the edge corresponding to halfedges `h` and `opposite(h,g)`, with the following invariant `halfedge(edge(h,g),g)==h`.
 */
template <typename HalfedgeGraph>
boost::graph_traits<HalfedgeGraph>::edge_descriptor
edge(boost::graph_traits<HalfedgeGraph>::halfedge_descriptor h, const HalfedgeGraph& g);

/*! \relates HalfedgeGraph
returns one of the halfedges corresponding to `e`.
 */
template <typename HalfedgeGraph>
boost::graph_traits<HalfedgeGraph>::halfedge_descriptor
halfedge(boost::graph_traits<HalfedgeGraph>::edge_descriptor f, const HalfedgeGraph& g);

/*! \relates HalfedgeGraph
returns a halfedge with target `v`.
 */
template <typename HalfedgeGraph>
boost::graph_traits<HalfedgeGraph>::halfedge_descriptor
halfedge(boost::graph_traits<HalfedgeGraph>::vertex_descriptor v, const HalfedgeGraph& g);


/*! \relates HalfedgeGraph
returns the halfedge with source `u` and target `v`. The Boolean is `true`, iff this halfedge exists.
 */
template <typename HalfedgeGraph>
std::pair<boost::graph_traits<HalfedgeGraph>::halfedge_descriptor,bool>
halfedge(boost::graph_traits<HalfedgeGraph>::vertex_descriptor u,
         boost::graph_traits<HalfedgeGraph>::vertex_descriptor v,
         const HalfedgeGraph& g);

/*! \relates HalfedgeGraph
returns the halfedge with source and target swapped.
 */
template <typename HalfedgeGraph>
boost::graph_traits<HalfedgeGraph>::halfedge_descriptor
opposite(boost::graph_traits<HalfedgeGraph>::halfedge_descriptor h, const HalfedgeGraph& g);

/*! \relates HalfedgeGraph
returns the source vertex of `h`.
 */
template <typename HalfedgeGraph>
boost::graph_traits<HalfedgeGraph>::vertex_descriptor
source(boost::graph_traits<HalfedgeGraph>::halfedge_descriptor h, const HalfedgeGraph& g);

/*! \relates HalfedgeGraph
returns the target vertex of `h`.
 */
template <typename HalfedgeGraph>
boost::graph_traits<HalfedgeGraph>::vertex_descriptor
target(boost::graph_traits<HalfedgeGraph>::halfedge_descriptor h, const HalfedgeGraph& g);

/*! \relates HalfedgeGraph
returns the next halfedge around its face.
 */
template <typename HalfedgeGraph>
boost::graph_traits<HalfedgeGraph>::halfedge_descriptor
next(boost::graph_traits<HalfedgeGraph>::halfedge_descriptor h, const HalfedgeGraph& g);

/*! \relates HalfedgeGraph
returns the previous halfedge around its face.
 */
template <typename HalfedgeGraph>
boost::graph_traits<HalfedgeGraph>::halfedge_descriptor
prev(boost::graph_traits<HalfedgeGraph>::halfedge_descriptor h, const HalfedgeGraph& g);

