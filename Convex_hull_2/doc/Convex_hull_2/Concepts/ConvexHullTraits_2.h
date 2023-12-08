/*!
\ingroup PkgConvexHull2Concepts
\cgalConcept

All convex hull and extreme point algorithms provided in \cgal are
parameterized with a traits class `Traits`, which defines the
primitives (objects and predicates) that the convex hull algorithms use.
`ConvexHullTraits_2` defines the complete set of primitives required in these
functions. The specific subset of these primitives required by each function
is specified with each function.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Convex_hull_constructive_traits_2<R>}
\cgalHasModels{CGAL::Convex_hull_traits_2<R>}
\cgalHasModels{CGAL::Convex_hull_traits_adapter_2<R>}
\cgalHasModels{CGAL::Projection_traits_xy_3<K>}
\cgalHasModels{CGAL::Projection_traits_yz_3<K>}
\cgalHasModels{CGAL::Projection_traits_xz_3<K>}
\cgalHasModelsEnd

*/

class ConvexHullTraits_2 {
public:

/// \name Types
/// @{

/*!
The point type on which the convex hull functions operate.
*/
typedef unspecified_type Point_2;

/*!
Binary predicate object type comparing `Point_2`s. Must provide
`bool operator()(Point_2 p, Point_2 q)` where `true`
is returned iff \f$ p ==_{xy} q\f$, false otherwise.

*/
typedef unspecified_type Equal_2;

/*!
Binary predicate object type comparing `Point_2`s
lexicographically. Must provide
`bool operator()(Point_2 p, Point_2 q)` where `true`
is returned iff \f$ p <_{xy} q\f$.
We have \f$ p<_{xy}q\f$, iff \f$ p_x < q_x\f$ or \f$ p_x = q_x\f$ and \f$ p_y < q_y\f$,
where \f$ p_x\f$ and \f$ p_y\f$ denote \f$ x\f$ and \f$ y\f$ coordinate of point \f$ p\f$,
respectively.

*/
typedef unspecified_type Less_xy_2;

/*!
Same as `Less_xy_2` with the roles of \f$ x\f$ and \f$ y\f$ interchanged.
*/
typedef unspecified_type Less_yx_2;

/*!
Predicate object type that must provide
`bool operator()(Point_2 p,Point_2 q,Point_2 r)`, which
returns `true` iff `r` lies to the left of the
oriented line through `p` and `q`.
*/
typedef unspecified_type Left_turn_2;

/*!
Predicate object type that must provide
`bool operator()(Point_2 p, Point_2 q, Point_2 r,Point_2 s)`,
which compares the signed distance of \f$ r\f$ and \f$ s\f$ to the directed line \f$ l_{pq}\f$
through \f$ p\f$ and \f$ q\f$.
It is used to compute the point right of a line with maximum unsigned distance to the line.
*/
typedef unspecified_type Compare_signed_distance_to_line_2;

/*!
Predicate object type that must provide
`bool operator()(Point_2 e, Point_2 p,Point_2 q)`,
where `true` is returned iff a tangent at \f$ e\f$ to the point set
\f$ \{e,p,q\}\f$ hits \f$ p\f$ before \f$ q\f$ when rotated counterclockwise around
\f$ e\f$.
Ties are broken such that the point with larger distance to \f$ e\f$
is smaller!
*/
typedef unspecified_type Less_rotate_ccw_2;

/*!
Predicate object type that must provide
`Orientation operator()(Point_2 e, Point_2 p,Point_2 q)`,
that returns `CGAL::LEFT_TURN`, if `r` lies to the left of the oriented line `l`
defined by `p` and `q`, returns `CGAL::RIGHT_TURN` if `r` lies to the right of `l`,
and returns `CGAL::COLLINEAR` if `r` lies on `l`.
*/
typedef unspecified_type Orientation_2;

/// @}

/// \name Creation
/// Only a copy constructor is required.
/// @{

/*!

*/
ConvexHullTraits_2(ConvexHullTraits_2& t);

/// @}

/// \name Operations
/// The following member functions to create instances of the above predicate
/// object types must exist.
/// @{

/*!

*/
Equal_2 equal_2_object();

/*!

*/
Less_xy_2 less_xy_2_object();

/*!

*/
Less_yx_2 less_yx_2_object();

/*!

*/
Compare_signed_distance_to_line_2 compare_signed_distance_to_line_2_object();

/*!

*/
Less_rotate_ccw_2 less_rotate_ccw_2_object( );

/*!

*/
Left_turn_2 left_turn_2_object();

/*!

*/
Orientation_2 orientation_2_object();

/// @}

}; /* end ConvexHullTraits_2 */
