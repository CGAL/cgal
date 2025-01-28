
/*!
\ingroup PkgAlphaShapes3Concepts
\cgalConcept

The concept `AlphaShapeCell_3` describes the requirements for the base cell of an alpha shape.

\cgalRefines{DelaunayTriangulationCellBase_3 if the underlying triangulation of the alpha shape is a Delaunay triangulation,
  RegularTriangulationCellBase_3 if the underlying triangulation of the alpha shape is a regular triangulation,
  Periodic_3TriangulationDSCellBase_3 if the underlying triangulation of the alpha shape is a periodic triangulation}

\cgalHasModelsBegin
\cgalHasModelsBare{`CGAL::Alpha_shape_cell_base_3` (templated with the appropriate triangulation cell base class)}
\cgalHasModelsEnd

\sa `CGAL::Alpha_status`

*/
class AlphaShapeCell_3 {
public:

/// \name Types
/// @{

/*!
A number type. Must be the same as the number type used
in the traits class of the triangulation underlying the alpha shape.
*/
typedef unspecified_type NT;

/*!
An iterator with value type `CGAL::Alpha_status<NT>`.
*/
typedef unspecified_type Alpha_status_iterator;

/// @}

/// \name Creation
/// @{

/*!
default constructor.
*/
AlphaShapeCell_3();

/*!
constructor setting the incident vertices.
*/
AlphaShapeCell_3(const Vertex_handle& v0, const Vertex_handle& v1, const Vertex_handle& v2, const Vertex_handle& v3);

/*!
constructor setting the incident vertices and the neighboring cells.
*/
AlphaShapeCell_3(const Vertex_handle& v0, const Vertex_handle& v1, const Vertex_handle& v2, const Vertex_handle& v3,
const Cell_handle& n0, const Cell_handle& n1, const Cell_handle& n2, const Cell_handle& n3);

/// @}

/// \name Access Functions
/// @{

/*!
Returns the alpha value of the cell.
*/
NT get_alpha();

/*!
Returns an iterator on the `CGAL::Alpha_status<NT>` of the facet
`i` of the cell.
*/
Alpha_status_iterator get_facet_status(int i);

/// @}

/// \name Modifiers
/// @{

/*!
Sets the critical value of the cell.
*/
void set_alpha(const NT & alpha);

/*!
Sets the iterator pointing to the `CGAL::Alpha_status<NT>`
of the facet `i` of the cell.
*/
void set_facet_status(int i, Alpha_status_iterator as);

/// @}

}; /* end AlphaShapeCell_3 */

