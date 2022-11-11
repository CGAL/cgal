
/*!
\ingroup PkgArrangementOnSurface2Concepts
\cgalConcept

A model of the `ArrangementVerticalRayShoot_2` concept can answer vertical ray-shooting
queries on an arrangement attached to it. Namely, given a
`Arrangement_2::Point_2` object, representing a point in the plane,
it returns the arrangement feature (edge or vertex) that lies
strictly above it (or below it). By "strictly" we mean that if the
query point lies on an arrangement edge (or on an arrangement vertex)
this edge will <I>not</I> be the query result, but the feature lying
above or below it. (An exception to this rule is the degenerate case
where the query point lies in the interior of a vertical edge.) Note
that it may happen that the query point lies above the upper envelope
(or below the lower envelope) of the arrangement, and the vertical ray
emanating from the query point goes to infinity without hitting any
arrangement feature on its way. In this case the unbounded face is
returned.

\cgalHeading{A Note on Backwards Compatibility}
The `ray_shoot_up` and `ray_shoot_down` member functions used
to return `CGAL::Object` up to \cgal version 4.2. Starting with
\cgal version 4.3 the return type is determined by a metafunction. To
preserve backwards compatibility `CGAL::Object` can be constructed
from the new return types implicitly, but switching to the new style
is recommended. To enable the old style without any overhead, the macro
`CGAL_ARR_POINT_LOCATION_VERSION` can be defined to 1 before any
\cgal header is included.

\cgalHasModel `CGAL::Arr_naive_point_location<Arrangement>`
\cgalHasModel `CGAL::Arr_walk_along_line_point_location<Arrangement>`
\cgalHasModel `CGAL::Arr_trapezoid_ric_point_location<Arrangement>`
\cgalHasModel `CGAL::Arr_landmarks_point_location<Arrangement,Generator>`

\sa `CGAL::Arr_naive_point_location<Arrangement>`
\sa `CGAL::Arr_walk_along_line_point_location<Arrangement>`
\sa `CGAL::Arr_trapezoid_ric_point_location<Arrangement>`
\sa `CGAL::Arr_landmarks_point_location<Arrangement,Generator>`
\sa `CGAL::Arr_point_location_result<Arrangement>`
\sa `CGAL_ARR_POINT_LOCATION_VERSION`
*/

class ArrangementVerticalRayShoot_2 {
public:

/// \name Types
/// @{

/*!
the associated arrangement type.
*/
typedef unspecified_type Arrangement_2;

/*!
equivalent to `Arrangement_2::Point_2`.
*/
typedef unspecified_type Point_2;

/// @}

/// \name Creation
/// @{

/*!
default constructor.
*/
ArrangementVerticalRayShoot_2();

/*!
constructs a ray-shooting object `rs` attached to the given
arrangement `arr`.
*/
ArrangementVerticalRayShoot_2 (const Arrangement_2& arr);

/// @}

/// \name Query Functions
/// @{

/*!
locates the arrangement feature that is first hit by an upward-directed
vertical ray emanating from the query point `q`,
and returns a handle for this feature. The function returns a
discriminated union container of the following bounded types:
<UL>
<LI>`Arrangement_2::Halfedge_const_handle`, in case the vertical
ray hits an arrangement edge;
<LI>`Arrangement_2::Vertex_const_handle`, in case the vertical
ray hits an arrangement vertex.
<LI>`Arrangement_2::Face_const_handle` for the unbounded arrangement
face, in case `q` lies above the upper envelope of the
arrangement.
</UL>
\pre `rs` is attached to a valid arrangement instance.
*/
Arr_point_location_result<Arrangement_2>::Type ray_shoot_up(const Point_2& q) const;

/*!
locates the arrangement feature that is first hit by a downward-directed
vertical ray emanating from the query point `q`,
and returns a handle for this feature. The function returns a
discriminated union container of the following bounded types:
<UL>
<LI>`Arrangement_2::Halfedge_const_handle`, in case the vertical
ray hits an arrangement edge;
<LI>`Arrangement_2::Vertex_const_handle`, in case the vertical
ray hits an arrangement vertex.
<LI>`Arrangement_2::Face_const_handle` for the unbounded arrangement
face, in case `q` lies below the lower envelope of the
arrangement.
</UL>
\pre `rs` is attached to a valid arrangement instance.
*/
Arr_point_location_result<Arrangement_2>::Type ray_shoot_down (const Point_2& q) const;

/// @}

/// \name Operations
/// @{

/*!
attaches `rs` to the given arrangement `arr`.
*/
void attach (const Arrangement_2& arr);

/*!
detaches `rs` from the arrangement it is currently attached to.
*/
void detach ();

/// @}

}; /* end ArrangementVerticalRayShoot_2 */

