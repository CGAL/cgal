
/*!
\ingroup PkgArrangement2Concepts
\cgalConcept

A model of the `ArrangementPointLocation_2` concept can answer point-location queries on
an arrangement attached to it. Namely, given a `Arrangement_2::Point_2`
object, representing a point in the plane, it returns the arrangement cell
containing it. In the general case, the query point is contained inside an
arrangement face, but in degenerate situations it may lie on an edge or
coincide with an arrangement vertex.

\cgalHeading{A note on Backwards compatibility}
The `locate` member function used to return `CGAL::Object` up to
\cgal version 4.2. Starting with \cgal version 4.3 the return type
is determined by a metafunction. To preserve backwards compatibility
`CGAL::Object` can be constructed from the new return types
implicitly, but switching to the new style is recommended. To enable
the old style without any overhead, the macro
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

class ArrangementPointLocation_2 {
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
ArrangementPointLocation_2(); 

/*!
constructs a point-location object `pl` attached to the given 
arrangement `arr`. 
*/ 
ArrangementPointLocation_2 (const Arrangement_2& arr); 

/// @} 

/// \name Query Functions 
/// @{

/*!
locates the arrangement cell that contains the query point `q` 
and returns a discriminated union container of the following bounded
types:

<UL> 
<LI>`Arrangement_2::Face_const_handle`, in case `q` is 
contained inside an arrangement face; 
<LI>`Arrangement_2::Halfedge_const_handle`, in case `q` lies 
on an arrangement edge; 
<LI>`Arrangement_2::Vertex_const_handle`, in case `q` coincides 
with an arrangement vertex. 
</UL> 
\pre `pl` is attached to a valid arrangement object. 
*/ 
Arr_point_location_result<Arrangement_2>::Type locate(const Point_2& q) const;

/// @} 

/// \name Operations 
/// @{

/*!
attaches `pl` to the given arrangement `arr`. 
*/ 
void attach (const Arrangement_2& arr); 

/*!
detaches `pl` from the arrangement it is currently attached to. 
*/ 
void detach (); 

/// @}

}; /* end ArrangementPointLocation_2 */

