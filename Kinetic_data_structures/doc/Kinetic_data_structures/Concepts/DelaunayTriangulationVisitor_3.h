namespace Kinetic {

/*!
\ingroup PkgKdsConcepts
\cgalConcept

This concept is for proxy objects which get notified when a kinetic Delaunay triangulation changes. 

\cgalHasModel `CGAL::Kinetic::Delaunay_triangulation_visitor_base_3`
\cgalHasModel `CGAL::Kinetic::Delaunay_triangulation_event_log_visitor_3` 

*/

class DelaunayTriangulationVisitor_3 {
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
void post_remove_vertex(Point_key, Cell_handle); 

/*!
The vertex is about to be inserted into the cell (although the `Cell_handle` might be `NULL`). 
*/ 
void pre_insert_vertex(Point_key, Cell_handle); 

/*!
The vertex was just inserted. 
*/ 
void post_insert_vertex(Vertex_handle); 

/*!
The trajectory of the point at the vertex changed. 
*/ 
void change_vertex(Vertex_handle); 

/*!
The cell was just created and initialized. 
*/ 
void create_cell(Cell_handle); 

/*!
The cell is about to be removed. 
*/ 
void destroy_cell(Cell_handle); 

/*!
The edge is about to be flipped. 
*/ 
void pre_edge_flip(Edge); 

/*!
The facet was just created with a flip. 
*/ 
void post_edge_flip(Facet); 

/*!
The facet is about to be flipped. 
*/ 
void pre_facet_flip(Facet); 

/*!
The edge was just created with a flip. 
*/ 
void post_facet_flip(Edge); 

/// @}

}; /* end DelaunayTriangulationVisitor_3 */

} /* end namespace Kinetic */
