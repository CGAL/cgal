namespace Kinetic {

/*!
\ingroup PkgKdsConcepts
\cgalConcept

This concept is for proxy objects which have functions called on them when a `Kinetic::Sort<Traits, Visitor>`. 

\cgalHasModel `CGAL::Kinetic::Sort_visitor_base`
\cgalHasModel `CGAL::Kinetic::Sort_event_log_visitor` 

*/

class SortVisitor {
public:

/// \name Operations 
/// @{

/*!
The vertex is about to be deleted. 
*/ 
void pre_remove_vertex(Vertex_handle); 

/*!
The vertex was just removed. 
*/ 
void post_remove_vertex(Point_key); 

/*!
The vertex is about to be inserted into the cell (although the `Cell_handle` might be `NULL`). 
*/ 
void pre_insert_vertex(Point_key); 

/*!
The vertex was just inserted. 
*/ 
void post_insert_vertex(Vertex_handle); 

/*!
Something changed at the vertex. 
*/ 
void change_vertex(Vertex_handle); 

/*!
The pair of vertices is about to be exchanged. 
*/ 
void pre_swap(Vertex_handle, Vertex_handle); 

/*!
The pair of vertices was just swapped. 
*/ 
void post_swap(Vertex_handle, Vertex_handle); 

/// @}

}; /* end Kinetic::SortVisitor */

} /* end namespae KineticConcepts */
