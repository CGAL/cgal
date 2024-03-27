
/*!
\ingroup PkgBooleanSetOperations2Concepts
\cgalConcept

A face record in a \em Dcel data structure used by the
`General_polygon_set_2` and `Polygon_set_2` template classes
to represent the underlying internal `Arrangement_2` data structure.

\cgalRefines{ArrangementDcelFace}

\cgalHasModelsBegin
\cgalHasModels{CGAL::Gps_face_base}
\cgalHasModelsEnd

\sa `ArrangementDcel`
\sa `ArrangementDcelVertex`
\sa `GeneralPolygonSetDcelHalfedge`

*/

class GeneralPolygonSetDcelFace {
public:

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

/*!
returns an id associated to the face using `set_id()`;
If no call to `set_id()` was done or if `reset_id()` was called,
the returned value is not specified and `id_not_set()` must return `true`.
*/
std::size_t id() const;

/*!
returns `true` if `set_id()` was not called or if `reset_id()` was called,
amd `false` otherwise.
*/
bool id_not_set() const;

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

/*!
resets the id associated to the face.
*/
void reset_id();

/*!
sets the id associated to the face.
i should be different from the largest possible value that `std::size_t`
can store.
*/
void set_id(std::size_t i);


/// @}

}; /* end GeneralPolygonSetDcelFace */
