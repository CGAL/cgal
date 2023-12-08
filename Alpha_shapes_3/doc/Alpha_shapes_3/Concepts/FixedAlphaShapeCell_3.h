
/*!
\ingroup PkgAlphaShapes3Concepts
\cgalConcept

The concept `FixedAlphaShapeCell_3` describes the requirements for the base cell of a alpha shape with a fixed value alpha.

\cgalRefines{DelaunayTriangulationCellBase_3 if the underlying triangulation of the alpha shape is a Delaunay triangulation,
  RegularTriangulationCellBase_3 if the underlying triangulation of the alpha shape is a regular triangulation,
  Periodic_3TriangulationDSCellBase_3 if the underlying triangulation of the alpha shape is a periodic triangulation}

\cgalHasModelsBegin
\cgalHasModelsBare{`CGAL::Fixed_alpha_shape_cell_base_3` (templated with the appropriate triangulation cell base class)}
\cgalHasModelsEnd
*/
class FixedAlphaShapeCell_3 {
public:

/// \name Creation
/// @{

/*!
default constructor.
*/
FixedAlphaShapeCell_3();

/*!
constructor setting the incident vertices.
*/
FixedAlphaShapeCell_3(const Vertex_handle& v0, const Vertex_handle& v1, const Vertex_handle& v2, const Vertex_handle& v3);

/*!
constructor setting the incident vertices and the neighboring cells.
*/
FixedAlphaShapeCell_3(const Vertex_handle& v0, const Vertex_handle& v1, const Vertex_handle& v2, const Vertex_handle& v3,
const Cell_handle& n0, const Cell_handle& n1, const Cell_handle& n2, const Cell_handle& n3);

/// @}

/// \name Access Functions
/// @{

/*!
Returns the classification of the cell.
*/
Classification_type get_classification_type();

/*!
Returns the classification of the i-th facet of the cell.
*/
Classification_type get_facet_classification_type(int i);

/// @}

/// \name Modifiers
/// @{

/*!
Sets the classification of the cell.
*/
void set_classification_type(Classification_type type);

/*!
Sets the classification of the i-th facet of the cell.
*/
void set_facet_classification_type(int i, Classification_type type);

/// @}

}; /* end FixedAlphaShapeCell_3 */

