
/*!
\ingroup PkgArrangementOnSurface2Concepts
\cgalConcept

A model for the `ArrangementOutputFormatter` concept supports a set of functions that enable
writing an arrangement to an output stream using a specific format.

\cgalHasModel `CGAL::Arr_text_formatter<Arrangement>`
\cgalHasModel `CGAL::Arr_face_extended_text_formatter<Arrangement>`
\cgalHasModel `CGAL::Arr_extended_dcel_text_formatter<Arrangement>`

*/

class ArrangementOutputFormatter {
public:

/// \name Types
/// @{

/*!
the type of arrangement to output.
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
typedef typename Arrangement_2::Vertex_const_handle Vertex_const_handle;

/*!

*/
typedef typename Arrangement_2::Halfedge_const_handle Halfedge_const_handle;

/*!

*/
typedef typename Arrangement_2::Face_const_handle Face_const_handle;

/// @}

/// \name Creation
/// @{

/*!
default constructor.
*/
Arr_out_formatter();

/*!
constructs a formatter that writes to `os`.
*/
Arr_out_formatter (std::ostream& os);

/*!
directs `outf` to write to `os`.
*/
void set_out (std::ostream& os);

/// @}

/// \name Access Functions
/// @{

/*!
returns the stream that `outf` writes to.
\pre `outf` is directed to a valid output stream.
*/
std::ostream& out ();

/// @}

/// \name Formatted Output Functions
/// @{

/*!
writes a message indicating the beginning of the arrangement.
*/
void write_arrangement_begin ();

/*!
writes a message indicating the end of the arrangement.
*/
void write_arrangement_end ();

/*!
writes a size value, preceded by a given label.
*/
void write_size (const char *label, Size size);

/*!
writes a message indicating the beginning of the vertex records.
*/
void write_vertices_begin();

/*!
writes a message indicating the end of the vertex records.
*/
void write_vertices_end();

/*!
writes a message indicating the beginning of the edge records.
*/
void write_edges_begin();

/*!
writes a message indicating the end of the edge records.
*/
void write_edges_end();

/*!
writes a message indicating the beginning of the face records.
*/
void write_faces_begin();

/*!
writes a message indicating the end of the face records.
*/
void write_faces_end();

/*!
writes a message indicating the beginning of a single vertex record.
*/
void write_vertex_begin();

/*!
writes a message indicating the end of a single vertex record.
*/
void write_vertex_end();

/*!
writes a vertex index.
*/
void write_vertex_index (std::size_t idx);

/*!
writes a point.
*/
void write_point (const Point_2& p);

/*!
writes the auxiliary data associated with the vertex.
*/
void write_vertex_data (Vertex_const_handle v);

/*!
writes a message indicating the beginning of a single edge record.
*/
void write_edge_begin();

/*!
writes a message indicating the end of a single edge record.
*/
void write_edge_end();

/*!
writes a halfedge index.
*/
void write_halfedge_index (std::size_t idx);

/*!
writes an \f$ x\f$-monotone curve.
*/
void write_x_monotone_curve (const X_monotone_curve_2& c);

/*!
writes the auxiliary data associated with the halfedge.
*/
void write_halfegde_data (Halfedge_const_handle he);

/*!
writes a message indicating the beginning of a single face record.
*/
void write_face_begin();

/*!
writes a message indicating the end of a single face record.
*/
void write_face_end();

/*!
writes a message indicating the beginning of the outer CCB of the current face.
*/
void write_outer_ccb_begin();

/*!
writes a message indicating the end of the outer CCB of the current face.
*/
void write_outer_ccb_end();

/*!
writes a message indicating the beginning of the container of holes inside the
current face.
*/
void write_holes_begin();

/*!
writes a message indicating the end of the container of holes inside the
current face.
*/
void write_holes_end();

/*!
writes a message indicating the beginning a connected component's boundary.
*/
void write_ccb_halfedges_begin();

/*!
writes a message indicating the end of a connected component's boundary.
*/
void write_ccb_halfedges_end();

/*!
writes a message indicating the beginning of the container of isolated vertices
inside the current face.
*/
void write_isolated_vertices_begin();

/*!
writes a message indicating the end of the container of isolated vertices inside
the current face.
*/
void write_isolated_vertices_end();

/*!
writes the auxiliary data associated with the face.
*/
void write_face_data (Face_const_handle f);

/// @}

}; /* end ArrangementOutputFormatter */

