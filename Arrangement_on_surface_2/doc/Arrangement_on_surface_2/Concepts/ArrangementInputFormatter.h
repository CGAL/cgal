
/*!
\ingroup PkgArrangement2Concepts
\cgalConcept

A model for the `ArrangementInputFormatter` concept supports a set of functions that enable 
reading an arrangement from an input stream using a specific format. 

\cgalHasModel `CGAL::Arr_text_formatter<Arrangement>` 
\cgalHasModel `CGAL::Arr_face_extended_text_formatter<Arrangement>` 
\cgalHasModel `CGAL::Arr_extended_dcel_text_formatter<Arrangement>` 

*/

class ArrangementInputFormatter {
public:

/// \name Types 
/// @{

/*!
the type of arrangement to input. 
*/ 
typedef unspecified_type Arrangement_2; 

/*!
the point type. 
*/ 
typedef typename Arrangement_2::Point_2 Point_2; 

/*!
the \f$ x\f$-monotone curve type. 
*/ 
typedef typename Arrangement_2::X_monotone_curve_2 X_monotone_curve_2; 

/*!

*/ 
typedef typename Arrangement_2::Size Size; 

/*!

*/ 
typedef typename Arrangement_2::Vertex_handle Vertex_handle; 

/*!

*/ 
typedef typename Arrangement_2::Halfedge_handle Halfedge_handle; 

/*!

*/ 
typedef typename Arrangement_2::Face_handle Face_handle; 

/// @} 

/// \name Creation 
/// @{

/*!
default constructor. 
*/ 
Arr_in_formatter(); 

/*!
constructs a formatter that reads from `is`. 
*/ 
Arr_in_formatter (std::istream& is); 

/*!
directs `inf` to read from `is`. 
*/ 
void set_in (std::istream& is); 

/// @} 

/// \name Access Functions 
/// @{

/*!
returns the stream that `inf` reads from. 
\pre `inf` is directed to a valid output stream. 
*/ 
std::istream& in (); 

/// @} 

/// \name Formatted Input Functions 
/// @{

/*!
reads a message indicating the beginning of the arrangement. 
*/ 
void read_arrangement_begin (); 

/*!
reads a message indicating the end of the arrangement. 
*/ 
void read_arrangement_end (); 

/*!
reads a size value, which is supposed to be preceded by the given label. 
*/ 
Size read_size (const char *label = NULL); 

/*!
reads a message indicating the beginning of the vertex records. 
*/ 
void read_vertices_begin(); 

/*!
reads a message indicating the end of the vertex records. 
*/ 
void read_vertices_end(); 

/*!
reads a message indicating the beginning of the edge records. 
*/ 
void read_edges_begin(); 

/*!
reads a message indicating the end of the edge records. 
*/ 
void read_edges_end(); 

/*!
reads a message indicating the beginning of the face records. 
*/ 
void read_faces_begin(); 

/*!
reads a message indicating the end of the face records. 
*/ 
void read_faces_end(); 

/*!
reads a message indicating the beginning of a single vertex record. 
*/ 
void read_vertex_begin(); 

/*!
reads a message indicating the end of a single vertex record. 
*/ 
void read_vertex_end(); 

/*!
reads and returns a vertex index. 
*/ 
std::size_t read_vertex_index (); 

/*!
reads a point. 
*/ 
void read_point (Point_2& p); 

/*!
reads an auxiliary vertex-data object and associates it with the vertex `v`. 
*/ 
void read_vertex_data (Vertex_handle v); 

/*!
reads a message indicating the beginning of a single edge record. 
*/ 
void read_edge_begin(); 

/*!
reads a message indicating the end of a single edge record. 
*/ 
void read_edge_end(); 

/*!
reads and returns halfedge index. 
*/ 
std::size_t read_halfedge_index (); 

/*!
reads an \f$ x\f$-monotone curve. 
*/ 
void read_x_monotone_curve (X_monotone_curve_2& c); 

/*!
reads an auxiliary halfedge-data object and associates it with the halfedge `he`. 
*/ 
void read_halfegde_data (Halfedge_handle he); 

/*!
reads a message indicating the beginning of a single face record. 
*/ 
void read_face_begin(); 

/*!
reads a message indicating the end of a single face record. 
*/ 
void read_face_end(); 

/*!
reads a message indicating the beginning of the outer CCB of the current face. 
*/ 
void read_outer_ccb_begin(); 

/*!
reads a message indicating the end of the outer CCB of the current face. 
*/ 
void read_outer_ccb_end(); 

/*!
reads a message indicating the beginning of the container of holes inside the 
current face. 
*/ 
void read_holes_begin(); 

/*!
reads a message indicating the end of the container of holes inside the 
current face. 
*/ 
void read_holes_end(); 

/*!
reads a message indicating the beginning of an inner CCB of the current face. 
*/ 
void read_inner_ccb_begin(); 

/*!
reads a message indicating the end of an inner CCB of the current face. 
*/ 
void read_inner_ccb_end(); 

/*!
reads a message indicating the beginning a connected component boundary. 
*/ 
void read_ccb_halfedges_begin(); 

/*!
reads a message indicating the end of a connected component boundary. 
*/ 
void read_ccb_halfedges_end(); 

/*!
reads a message indicating the beginning of the container of isolated vertices 
inside the current face. 
*/ 
void read_isolated_vertices_begin(); 

/*!
reads a message indicating the end of the container of isolated vertices inside 
the current face. 
*/ 
void read_isolated_vertices_end(); 

/*!
reads an auxiliary face-data object and associates it with the face `f`. 
*/ 
void read_face_data (Face_handle f); 

/// @}

}; /* end ArrangementInputFormatter */

