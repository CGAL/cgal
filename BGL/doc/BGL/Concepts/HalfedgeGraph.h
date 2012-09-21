
/*!
\ingroup PkgBGLConcepts
\cgalconcept

The concept `HalfedgeGraph` describes the requirements for a graph that is 
structurally equivalent to a polyhedral surface represented by a 
halfedge data structure, and it provides an interface for efficient 
access to the opposite edge of an edge, and to the successor and 
predecessor of an edge in the iterator range of the incoming edges of 
a vertex. Each vertex has a geometric position in space. As in a 
halfedge data structure we define the face adjacent to a halfedge to be 
to the <I>left</I> of the halfedge. 

### Requirements ###

For each <I>directed edge</I> \f$ e=(v,w)\f$ its opposite edge \f$ e'=(w,v)\f$ 
must be part of the graph. 

The incoming edges of a vertex \f$ v\f$ have a fixed order, that is all 
calls of `in_edges(v,g)` must return the same iterator range, 
modulo a cyclic permutation. The order must be <I>clockwise</I>. 

As the `HalfedgeGraph` is equivalent to a polyhedral surface there must exist an embedding 
for the vertices and edges such that the ordered edges do not intersect. 

\refines <A HREF="http://www.boost.org/libs/graph/doc/BidirectionalGraph.html">BidirectionalGraph</A> 
\refines <A HREF="http://www.boost.org/libs/graph/doc/PropertyGraph.html">PropertyGraph</A> 

A model of `HalfedgeGraph` must have the <I>interior properties</I>
`edge_is_border` attached to its edges, and it must have
`vertex_is_border` and `vertex_point` attached to its vertices.

### Associated Types ###

Because (directed) edges must come in pairs, there is the
additional notion of an <I>undirected edge</I>
\footnote{The directed edges are not called `halfedges` (as in a `HalfedgeDS`) because from the point of view of this graph, being a refinement of a sc{Bgl} graph, each directed edge is an edge in itself. In other words, the unqualified term edge refers to one and only one directed edge and not to a pair.} 
for a pair of opposite directed edges. The number of undirected
edges is exactly half the number of directed edges. Note that the
notion of directed and undirected edges does not imply the
existence of two different types. The type `edge_descriptor` is
used for both. An undirected edge must be implicitly handled, and
there is no requirement on which of the directed edges of the
undirected edge must be used to represent it.

\code
// The type of the geometric location of a vertex. 
typedef Hidden_type halfedge_graph_traits<HalfedgeGraph>::Point; 

// An iterator that iterates over one and only one of the directed edges 
// in each pair of opposite directed edges. The value type of the iterator 
// is `boost::graph_traits<HalfedgeGraph>::edge_descriptor`. 
typedef Hidden_type halfedge_graph_traits<HalfedgeGraph>::undirected_edge_iterator; 
\endcode

### Valid Expressions ###

Following the \sc{Bgl} design, the following graph operations are defined as free 
rather than member functions. 

\hasModel `CGAL::Polyhedron_3<Traits>`
*/
class HalfedgeGraph {

};

/*! 
Returns the undirected edges of `g`. 

\relates HalfedgeGraph 
*/ 
template<class Graph> 
std::pair<typename halfedge_graph_traits<HalfedgeGraph>::undirected_edge_iterator, 
          typename halfedge_graph_traits<HalfedgeGraph>::undirected_edge_iterator> 
undirected_edges(const Graph& g); 

/*! 
Returns the opposite edge of `e`. 

An edge \f$ e=(v,w)\f$ is said to be the <I>opposite edge</I> of edge
\f$ e'=(w,v)\f$.

\relates HalfedgeGraph 
*/ 
template<class Graph> 
typename boost::graph_traits<Graph const>::edge_descriptor 
opposite_edge(typename boost::graph_traits<Graph const>::edge_descriptor e, Graph const& g ); 

/*! 
Returns the clockwise neighbor of `e`. 

An edge \f$ e'=(v,w)\f$ is called the <I>clockwise neighbor</I> of
edge \f$ e=(u,w)\f$, and \f$ e\f$ the <I>counterclockwise neighbor</I>
of \f$ e'\f$, iff there exist two iterators \f$ it\f$ and \f$ it'\f$
in the iterator range `in_edges(w,g)` such that `**it == e` and `**it'
== e'`, and `it' == it++` or `it` is the last and `it'` the first
iterator of the iterator range.

\relates HalfedgeGraph 
*/ 
template<class Graph> 
typename boost::graph_traits<Graph const>::edge_descriptor 
next_edge_cw(typename boost::graph_traits<Graph const>::edge_descriptor e, Graph const& g ); 

/*! 
Returns the counterclockwise neighbor of `e`. 
\relates HalfedgeGraph 
*/ 
template<class Graph> 
typename boost::graph_traits<Graph const>::edge_descriptor 
next_edge_ccw(typename boost::graph_traits<Graph const>::edge_descriptor e, Graph const& g ); 

/*! 
Returns the successor of `e`. 

An edge \f$ e'=(v,w)\f$ is called the <I>successor</I> of edge \f$
e=(u,v)\f$, and \f$ e\f$ the <I>predecessor</I> of \f$ e'\f$, iff \f$
e'\f$ is the clockwise neighbor of the opposite edge of \f$ e\f$.

\relates HalfedgeGraph 
*/ 
template<class Graph> 
typename boost::graph_traits<Graph const>::edge_descriptor 
next_edge(typename boost::graph_traits<Graph const>::edge_descriptor e, Graph const& g ); 

/*! 
Returns the predecessor of `e`. 

An edge \f$ e'=(v,w)\f$ is called the <I>successor</I> of edge \f$
e=(u,v)\f$, and \f$ e\f$ the <I>predecessor</I> of \f$ e'\f$, iff \f$
e'\f$ is the clockwise neighbor of the opposite edge of \f$ e\f$.

\relates HalfedgeGraph 
*/ 
template<class Graph> 
typename boost::graph_traits<Graph const>::edge_descriptor 
prev_edge(typename boost::graph_traits<Graph const>::edge_descriptor e, Graph const& g );
