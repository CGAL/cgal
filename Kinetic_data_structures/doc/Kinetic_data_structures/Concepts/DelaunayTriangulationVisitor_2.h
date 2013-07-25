namespace Kinetic {

/*!
\ingroup PkgKdsConcepts
\cgalConcept

This concept is for proxy objects which get notified when a kinetic Delaunay triangulation changes. 

\cgalHasModel `CGAL::Kinetic::Delaunay_triangulation_visitor_base_2`
\cgalHasModel `CGAL::Kinetic::Delaunay_triangulation_recent_edges_visitor_2<Triangulation>`
\cgalHasModel `CGAL::Kinetic::Delaunay_triangulation_event_log_visitor_2 `

*/

class DelaunayTriangulationVisitor_2 {
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
The trajectory of the point at the vertex changed. 
*/ 
void change_vertex(Vertex_handle); 

/*!
The face has just been made. 
*/ 
void create_face(Face_handle f); 

/*!
The face is about to be destroyed. 
*/ 
void destroy_face(Face_handle f); 

/*!
The edge is about to be flipped. 
*/ 
void pre_flip(Edge); 

/*!
The edge was just created with a flip. 
*/ 
void post_flip(Edge); 

/// @}

}; /* end DelaunayTriangulationVisitor_2 */

} /* end namespace Kinetic */
