namespace boost {

/*!
\ingroup PkgBGLTraits

The class `graph_traits` is a partial specialization of 
<A HREF="http://www.boost.org/libs/graph/doc/graph_traits.html">`boost::graph_traits`</A> 
for the 2D triangulation classes.

The triangulations of \cgal are all models of the concepts
`BidirectionalGraph` and `VertexAndEdgeListGraph` of the Boost Graph
Library \cite cgal:sll-bgl-02.


The mapping between vertices and edges of the triangulation and the 
graph is rather straightforward, but there are some subtleties. The 
value type of the \sc{Bgl} iterators is the vertex or edge descriptor, 
whereas in \cgal all iterators and circulators are also handles and 
hence have as value type Vertex or Edge. 

The graph traits class for triangulations does not distinguish between 
finite and infinite vertices and edges. As the edge weight computed 
with the default property map of \sc{Bgl} algorithms (obtained with 
`boost::get(t, boost::edge_weight)`) is the length of the edge, 
the edge weight is not well defined for infinite edges. For algorithms 
that make use of the edge weight the user must therefore 
define a <A HREF="http://www.boost.org/libs/graph/doc/filtered_graph.html">`boost::filtered_graph`</A> or pass a property map to the 
algorithm that returns "infinity" for infinite edges. 

Note also that when you derive from the class `CGAL::Triangulation_2` 
you must upcast the object in order to use this partial specialization. 

For the user convenience, \cgal provides the partial specializations 
for all 2D triangulation classes. 

*/
template< typename GT, typename TDS> >
class graph_traits< CGAL::Triangulation_2<GT, TDS> > {
public:

/// \name Types 
/// @{

/*! 
The vertex descriptor. 
*/ 
typedef Triangulation::Vertex_handle vertex_descriptor; 

/*! 
The edge descriptor. It is constructible from and convertible to `Triangulation::Edge`. 
The edge descriptor is not a simple typedef, but a proper class, 
because in an undirected graph 
the edges `(u,v)` and `(v,u)` must be equal. This is not the case 
for the Edge type of the triangulation. 
*/ 
typedef Hidden_type edge_descriptor; 

/*! 
The vertex iterator type. Its value type is `vertex_descriptor`. 
*/ 
typedef Hidden_type vertex_iterator; 

/*! 
The edge iterator type, Its value type is `edge_descriptor`. 
*/ 
typedef Hidden_type edge_iterator; 

/*! 
An iterator for the outgoing edges incident to a vertex. 
Its value type is `edge_descriptor`. 
*/ 
typedef Hidden_type out_edge_iterator; 

/*! 
An iterator for the incoming edges incident to a vertex. 
Its value type is `edge_descriptor`. 
*/ 
typedef Hidden_type in_edge_iterator; 

/*! 
An iterator for the vertices adjacent to a vertex. 
Its value type is `vertex_descriptor`. 
*/ 
typedef Hidden_type adjacency_iterator; 

/*! 

*/ 
typedef boost::undirected_tag directed_category; 

/*! 

*/ 
typedef boost::disallow_parallel_edge_tag edge_parallel_category; 

/// @}

}; /* end graph_traits */
} /* end namespace boost */
