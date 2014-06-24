
namespace CGAL {

/*!
\ingroup PkgBGLTraits

\deprecated This class is deprecated since \cgal 4.5, the class
`boost::graph_traits` should be used instead.
 
The class `halfedge_graph_traits` is a partial specialization of `CGAL::halfedge_graph_traits` 
for `Polyhedron_3`. It provides the types associated 
to the `HalfedgeGraph` concept. 

*/
template< typename Traits> >
class halfedge_graph_traits< Polyhedron_3<Traits> > {
public:

/// \name Types 
/// @{

/*!
An edge iterator that iterates 
over one of the two opposite edges forming an undirected edge. 

The value type is `CGAL::Polyhedron_3::Halfedge_const_handle`. 
*/ 
typedef unspecified_type undirected_edge_iterator; 

/*!
The point type of the vertex. 
*/ 
typename CGAL::Polyhedron_3<Traits>::Point_3 Point ; 

/// @}

}; /* end halfedge_graph_traits */
} /* end namespace CGAL */
