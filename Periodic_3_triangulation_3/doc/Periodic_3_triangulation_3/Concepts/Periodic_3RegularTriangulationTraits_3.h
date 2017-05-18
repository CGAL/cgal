
/*!
\ingroup PkgPeriodic3Triangulation3Concepts
\cgalConcept

The concept `Periodic_3RegularTriangulationTraits_3` is the first template parameter
of the class `CGAL::Periodic_3_regular_triangulation_3`.
It refines the concept
`RegularTriangulationTraits_3` from the \cgal 3D Triangulations.
It redefines the geometric objects, predicates and constructions to
work with point-offset pairs. In most cases the offsets will be
(0,0,0) and the predicates from `RegularTriangulationTraits_3`
can be used directly. For efficiency reasons we maintain for each
functor the version without offsets.

\cgalRefines Periodic_3TriangulationTraits_3
\cgalRefines RegularTriangulationTraits_3

\cgalHasModel CGAL::Periodic_3_regular_triangulation_traits_3

In addition to the requirements described for the traits class
RegularTriangulationTraits_3, the geometric traits class of a
periodic regular triangulation must fulfill the following
requirements.

\note The optional types must be provided in any case, however they
can be replaced by dummy types if the respective functions are not
used.
*/
class Periodic_3RegularTriangulationTraits_3 {
public:

/// \name
/// @{

/*!
A predicate object that must provide the function operators:

`Oriented_side operator()(Weighted_point_3 p, Weighted_point_3 q, Weighted_point_3 r, Weighted_point_3 s, Weighted_point_3 t)`,

which determines the position of `t` with respect to the power sphere of `p, q, r, s`.

`Oriented_side operator()(Weighted_point_3 p, Weighted_point_3 q, Weighted_point_3 r, Weighted_point_3 s, Weighted_point_3 t,
Periodic_3_offset_3 o_p, Periodic_3_offset_3 o_q, Periodic_3_offset_3 o_r, Periodic_3_offset_3 o_s, Periodic_3_offset_3 o_t)`,

which is the same for the point-offset pair `(t,o_t)` with respect to the power sphere of the point-offset pairs
`(p,o_p), (q,o_q), (r,o_r), (s,o_s)`.
\pre `p`, `q`, `r`, `s`, `t` lie inside the domain and `p, q, r, s` are not coplanar.

<HR WIDTH=50%>

When vertex removal is used, the predicate must in addition provide the following function operators:

`Oriented_side operator()( Weighted_point_3 p, Weighted_point_3 q, Weighted_point_3 r, Weighted_point_3 t)`,

which has a definition similar to the previous method, for coplanar points,
with the power circle of `p,q,r`.

`Oriented_side operator()( Weighted_point_3 p, Weighted_point_3 q, Weighted_point_3 r, Weighted_point_3 t,
Periodic_3_offset_3 o_p, Periodic_3_offset_3 o_q, Periodic_3_offset_3 o_r, Periodic_3_offset_3 o_t)`,

which is the same for point-offset pairs.

\pre `p`, `q`, `r`, `t` lie inside the domain, `p, q, r` are not collinear, and `p, q, r, t` are coplanar.


`Oriented_side operator()( Weighted_point_3 p, Weighted_point_3 q, Weighted_point_3 t)`,

which is the same for collinear points, and the power segment of `p` and `q`,

`Oriented_side operator()( Weighted_point_3 p, Weighted_point_3 q, Weighted_point_3 t,
Periodic_3_offset_3 o_p, Periodic_3_offset_3 o_q, Periodic_3_offset_3 o_t)`,

which is the same for point-offset pairs.

\pre `p`, `q`, `t` lie inside the domain, `p` and `q` have different Bare_points, and `p, q, t` are collinear.


`Oriented_side operator()( Weighted_point_3 p, Weighted_point_3 q)`,

which is the same for equal points, that is when `p` and `q`
have equal coordinates, then it returns the comparison of the weights.

`Oriented_side operator()( Weighted_point_3 p, Weighted_point_3 q,
Periodic_3_offset_3 o_p, Periodic_3_offset_3 o_q)`,

which is the same for point-offset pairs.

\pre `p` and `q` lie inside the domain and have equal Bare_points.

*/
typedef unspecified_type Power_side_of_oriented_power_sphere_3;

/// @}

/// \name
/// @{

/*!
A predicate object that must provide the function operators:

`Orientation operator()(Weighted_point_3 p, Weighted_point_3 q, Weighted_point_3 r, Weighted_point_3 s, FT w)`,

which compares the weight of the smallest sphere orthogonal to the input weighted
points with the input weight `w` and returns a `SMALLER`, `EQUAL`, or `LARGER`.

\pre `p`, `q`, `r`, and `s` lie inside the domain.

*/
typedef unspecified_type Compare_weighted_squared_radius_3;

/// @}

/// \name
/// @{

/*!
A predicate object, model of `ComparePowerDistance_3`, that must provide
the function operator

`Comparison_result operator()(Point_3 p, Weighted_point_3 q, Weighted_point_3 r)`,

which compares the power distance between `p` and `q` to the power distance
between `p` and `r` and

`Comparison_result operator()(Point_3 p, Weighted_point_3 q, Weighted_point_3 r,
Periodic_3_offset_3 o_p, Periodic_3_offset_3 o_q, Periodic_3_offset_3 o_r)`,

which is the same for point-offset pairs.

\note This predicate is required if a call to `nearest_power_vertex()` or
`nearest_power_vertex_in_cell()` is issued.*/
typedef unspecified_type Compare_power_distance_3;

/// @}

/// \name
/// When vertex removal is used, the traits class must in addition provide the following predicate object
/// @{

/*!
A predicate object that must provide the function operators:

`Orientation operator()(Weighted_point_3 p, Weighted_point_3 q, Weighted_point_3 r)`,

which returns `COLLINEAR`, if the points are collinear; otherwise
it must return a consistent orientation for any three points chosen in
a same plane and

`Orientation operator()(Weighted_point_3 p, Weighted_point_3 q, Weighted_point_3 r,
Periodic_3_offset_3 o_p, Periodic_3_offset_3 o_q, Periodic_3_offset_3 o_r)`,

which is the same for point-offset pairs.
\pre `p`, `q`, `r` lie inside the domain.
*/
typedef unspecified_type Coplanar_orientation_3;

/// @}

/// \name
/// When the dual operations are used, the traits
/// class must in addition provide the following constructor object:
/// @{

/*!
A constructor object that must provide the function operators:

`Weighted_point_3 operator()(Weighted_point_3 p, Weighted_point_3 q, Weighted_point_3 r, Weighted_point_3 s)`,

which constructs the weighted circumcenter of four points and

`Weighted_point_3 operator()(Weighted_point_3 p, Weighted_point_3 q, Weighted_point_3 r, Weighted_point_3 s,
Periodic_3_offset_3 o_p, Periodic_3_offset_3 o_q, Periodic_3_offset_3 o_r, Periodic_3_offset_3 o_s)`,

which constructs the weighted circumcenter of four point-offset pairs.
\pre `p`, `q`, `r` and `s` as well as `(p,o_p)`, `(q,o_q)`, `(r,o_r)` and `(s,o_s)` must be non coplanar. `p`, `q`, `r`, `s` lie inside the domain.
*/
typedef unspecified_type Construct_weighted_circumcenter_3;

/// @}

/// \name Operations
/// The following functions give access to the predicate and construction objects:
/// @{

/*!

*/
Power_side_of_oriented_power_sphere_3 power_side_of_oriented_power_sphere_3_object();

Compare_weighted_squared_radius_3 compare_weighted_squared_radius_3_object();

/// @}

/// \name
/// The following function must be provided if vertex removal is
/// used; otherwise dummy functions can be provided.
/// @{

/*!

*/
Coplanar_orientation_3 coplanar_3_orientation_3_object();

/// @}

/// \name
/// The following function must be provided only if the methods of
/// `Periodic_3_regular_triangulation_3` returning elements of the
/// Voronoi diagram are used; otherwise a dummy function can be
/// provided.
/// @{

/*!

*/
Construct_weighted_circumcenter_3 construct_weighted_circumcenter_3_object();

/// @}

}; /* end Periodic_3RegularTriangulationTraits_3 */

