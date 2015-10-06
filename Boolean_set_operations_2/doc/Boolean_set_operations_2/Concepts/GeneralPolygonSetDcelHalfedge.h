/*!
\ingroup PkgBooleanSetOperations2Concepts
\cgalConcept

A halfedge record in a \em Dcel data structure used
by the `General_polygon_set_2` and `Polygon_set_2` template classes
to represent the underlying internal `Arrangement_2` data structure.

\cgalRefines `ArrangementDcelHalfedge`

\cgalHasModel `CGAL::Gps_face_halfedge`

\sa `ArrangementDcel`
\sa `ArrangementDcelVertex`
\sa `GeneralPolygonSetDcelFace`

*/

class GeneralPolygonSetDcelHalfedge {
public:

/// \name Access Functions
/// @{
/*!
returns an integer passed to `set_flag()`,
if the latter function was not called the output
is not specified.
*/
int flag() const;


/// @}

/// \name Modifiers
/// @{
/*!
set an integer flag for this halfedge.
*/
void set_flag(int i);
/// @}

}; /* end GeneralPolygonSetDcelHalfedge */
