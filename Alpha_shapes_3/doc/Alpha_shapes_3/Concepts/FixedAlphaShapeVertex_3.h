
/*!
\ingroup PkgAlphaShapes3Concepts
\cgalConcept

The concept `FixedAlphaShapeVertex_3` describes the requirements for the base vertex of a alpha shape with a fixed value alpha.

\cgalRefines `TriangulationVertexBase_3`, if the underlying triangulation of the alpha shape is a Delaunay triangulation.
\cgalRefines `RegularTriangulationVertexBase_3`, if the underlying triangulation of the alpha shape is a regular triangulation.
\cgalRefines `Periodic_3TriangulationDSVertexBase_3`, if the underlying triangulation of the alpha shape is a periodic triangulation.

\cgalHasModel `CGAL::Fixed_alpha_shape_vertex_base_3` (templated with the appropriate triangulation vertex base class).
*/

class FixedAlphaShapeVertex_3 {
public:

/// \name Types
/// @{

/*!
Must be the same as the point type provided by
the geometric traits class of the triangulation.
*/
typedef unspecified_type Point;

/// @}

/// \name Creation
/// @{

/*!
default constructor.
*/
FixedAlphaShapeVertex_3();

/*!
constructor setting the point.
*/
FixedAlphaShapeVertex_3(Point p);

/*!
constructor setting the point and an incident cell.
*/
FixedAlphaShapeVertex_3(Point p, const Cell_handle& c);

/// @}

/// \name Access Functions
/// @{

/*!
Returns a boolean indicating whether the point is on the convex hull of the point of the triangulation.
*/
bool is_on_chull();

/*!
Returns the classification of the vertex.
*/
Classification_type get_classification_type();

/// @}

/// \name Modifiers
/// @{

/*!
Sets the classification of the vertex.
*/
void set_classification_type(Classification_type type);

/*!
Sets whether the vertex is on the convex hull.
*/
void is_on_chull(bool b);

/// @}

}; /* end FixedAlphaShapeVertex_3 */

