/*!
\ingroup PkgBooleanSetOperations2Concepts
\cgalConcept

A halfedge record in a <span class="textsc">Dcel</span> data structure used
by the `General_polygon_set_2` and `Polygon_set_2` class-templates
to represent the underlying internal `Arrangement_2` data structure.

\cgalRefines `ArrangementDcelHalfedge`

\sa `ArrangementDcel`
\sa `ArrangementDcelVertex`
\sa `GeneralPolygonSetDcelFace`

*/

class GeneralPolygonSetDcelHalfedge {
public:

/// \name Creation
/// @{

/*!
default constructor.
*/
Gps_dcel_halfedge();

/*!
assigns `f` with the contents of the `other` face.
*/
void assign(const Self& other);

/// @}

/// \name Access Functions
/// @{

/// @}

/// \name Modifiers
/// @{

/// @}

}; /* end GeneralPolygonSetDcelHalfedge */
