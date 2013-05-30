namespace boost {

/*!
\ingroup PkgBGLTraits

The class `graph_traits` is a partial specialization of 
<A HREF="http://www.boost.org/libs/graph/doc/graph_traits.html">`boost::graph_traits`</A> 
for the class `CGAL::Arrangement_2`. It provides the types associated 
to the 
<A HREF="http://www.boost.org/libs/graph/doc/graph_concepts.html">graph</A> concepts 
<A HREF="http://www.boost.org/libs/graph/doc/BidirectionalGraph.html">`BidirectionalGraph`</A> and 
<A HREF="http://www.boost.org/libs/graph/doc/EdgeAndVertexListGraph.html">`EdgeAndVertexListGraph`</A>. 

The const specialization, `boost::graph_traits< CGAL::Arrangement_2<Traits,Dcel> const>` 
is also defined, using the constant handles in the arrangement. 

*/
template< typename T, typename DC> >
class graph_traits< CGAL::Arrangement_2<T, DC> > {
public:

/// \name Types 
/// @{

/*! 
The vertex descriptor. 
*/ 
typename CGAL::Arrangement_2::Vertex_handle vertex_descriptor; 

/*! 
The edge descriptor. 
*/ 
typename CGAL::Arrangement_2::Halfedge_handle edge_descriptor; 

/*! 
An iterator corresponding to 
`CGAL::Arrangement_2::Vertex_iterator`, 
with the difference that its value type is a vertex descriptor and not 
`CGAL::Arrangement_2::Vertex`. 
*/ 
typedef Hidden_type vertex_iterator; 

/*! 
An iterator corresponding to 
`CGAL::Arrangement_2::Halfedge_iterator` 
with the difference that its value type is an edge descriptor and not 
`CGAL::Arrangement_2::Halfedge`. 
*/ 
typedef Hidden_type edge_iterator; 

/*! 
An edge iterator which only iterates over 
the incoming edges around a vertex. It corresponds to a 
`CGAL::Arrangement_2::Halfedge_around_vertex_circulator` 
with the difference that its value type is an edge descriptor and not 
`CGAL::Arrangement_2::Halfedge`. 
*/ 
typedef Hidden_type in_edge_iterator; 

/*! 
An edge iterator which only iterates over 
the outgoing halfedges around a vertex. It corresponds to a 
`CGAL::Arrangement_2::Halfedge_around_vertex_circulator` 
with the difference that its value type is an edge descriptor and not 
`CGAL::Arrangement_2::Halfedge`. 
*/ 
typedef Hidden_type out_edge_iterator; 

/*! 
Indicates that this graph does support multiedges. 
*/ 
boost::disallow_parallel_edge_tag edge_parallel_category; 

/*! 
Indicates that this graph is bidirectional. 
*/ 
boost::bidirectional_graph_tag traversal_category; 

/*! 
The size type of the vertex list. 
*/ 
typename Arrangement_2::size_type vertices_size_type; 

/*! 
The size type of the edge list. 
*/ 
typename Arrangement_2::size_type edges_size_type; 

/*! 
The size type of the adjacency list. 
*/ 
typename Arrangement_2::size_type degree_size_type; 

/// @}

}; /* end graph_traits */
} /* end namespace boost */
