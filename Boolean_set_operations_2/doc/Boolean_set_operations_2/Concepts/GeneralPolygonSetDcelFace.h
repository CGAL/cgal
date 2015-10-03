
/*!
\ingroup PkgBooleanSetOperations2Concepts
\cgalConcept

A face record in a <span class="textsc">Dcel</span> data structure used by the
`General_polygon_set_2` and `Polygon_set_2` class-templates
to represent the underlying internal `Arrangement_2` data structure.

\cgalRefines `ArrangementDcelFace`

\sa `ArrangementDcel`
\sa `ArrangementDcelVertex`
\sa `GeneralPolygonSetDcelHalfedge`

*/

class GeneralPolygonSetDcelFace {
public:

/// \name Creation
/// @{

/*!
default constructor.
*/
Gps_dcel_face();

/*!
assigns `f` with the contents of the `other` face.
*/
void assign (const Self& other);

/// @}

/// \name Access Functions
/// @{

/*!
returns `true` if the face is contained in the general-polygon set,
and `false` otherwise.
*/
bool contained() const;

/*!
returns `true` if the face has been visited, and `false` otherwise.
This is used internally by the some of the operations of the
`General_polygon_set_2` class that traverse the arrangement faces.
*/
bool visited();

/// @}

/// \name Modifiers
/// @{

/*!
marks the face as contained (if `flag` is `true`), or as a hole
(if it is `false`).
*/
void set_contained(bool flag);

/*!
marks the face as visited (if `flag` is `true`), or as not visited
(if it is `false`). This is used internally by the some of the
operations of the `General_polygon_set_2` class that traverse the
arrangement faces.
*/
void set_visited(bool flag);

/// @}

}; /* end GeneralPolygonSetDcelFace */
