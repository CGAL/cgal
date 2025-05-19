
/*!
\ingroup PkgPeriodic3Triangulation3Concepts
\cgalConcept

The concept `Periodic_3DelaunayTriangulationTraits_3` is the first template parameter
of the class `CGAL::Periodic_3_Delaunay_triangulation_3`.
It refines the concepts `Periodic_3TriangulationTraits_3` and
`DelaunayTriangulationTraits_3`.
It redefines the geometric objects, predicates and constructions to
work with point-offset pairs. In most cases the offsets will be
(0,0,0) and the predicates from `DelaunayTriangulationTraits_3`
can be used directly. For efficiency reasons we maintain for each
functor the version without offsets.

\cgalRefines{Periodic_3TriangulationTraits_3,DelaunayTriangulationTraits_3}

\cgalHasModelsBegin
\cgalHasModels{CGAL::Periodic_3_Delaunay_triangulation_traits_3}
\cgalHasModelsEnd

In addition to the requirements described by the concepts
`Periodic_3TriangulationTraits_3` and `DelaunayTriangulationTraits_3`,
the geometric traits class of a Periodic Delaunay triangulation must fulfill
the following requirements.

\note The optional types must be provided in any case, however they
can be replaced by dummy types if the respective functions are not
used.
*/
class Periodic_3DelaunayTriangulationTraits_3 {
public:

/// \name
/// @{

/*!
A predicate object that must provide the function operator

`Oriented_side operator()(Point_3 p, Point_3 q, Point_3 r, Point_3 s, Point_3 t, Periodic_3_offset_3 o_p, Periodic_3_offset_3 o_q, Periodic_3_offset_3 o_r, Periodic_3_offset_3 o_s, Periodic_3_offset_3 o_t)`,

which determines on which side of the oriented sphere circumscribing
`(p,o_p), (q,o_q), (r,o_r), (s,o_s)` the point-offset pair
`(t,o_t)` lies.
\pre `p`, `q`, `r`, `s`, `t` lie inside the domain.
*/
typedef unspecified_type Side_of_oriented_sphere_3;

/*!
A predicate object that must provide the function operator

`Comparison_result operator()(Point_3 p, Point_3 q, Point_3 r, Periodic_3_offset_3 o_p, Periodic_3_offset_3 o_q, Periodic_3_offset_3 o_r)`,

which compares the distance between `(p,o_p)` and `(q,o_q)` to
the distance between `(p,o_p)` and `(r,o_r)`.
\pre `p`, `q`, `r` lie inside the domain.
*/
typedef unspecified_type Compare_distance_3;

/// @}

/// \name
/// When vertex removal is used, the traits class must in addition provide the following predicate objects
/// @{

/*!
A predicate object that must provide the function operator

`Orientation operator()(Point_3 p, Point_3 q, Point_3 r Periodic_3_offset_3 o_p, Periodic_3_offset_3 o_q, Periodic_3_offset_3 o_r)`,

which returns `COLLINEAR`, if the point-offset pairs are
collinear; otherwise it must return a consistent orientation for any
three point-offset pairs chosen in a same plane.
\pre `p`, `q`, `r` lie inside the domain.
*/
typedef unspecified_type Coplanar_orientation_3;

/*!
A predicate object that must provide the function operator

`Bounded_side operator()(Point_3 p, Point_3 q, Point_3 r, Point_3 s, Periodic_3_offset_3 o_p, Periodic_3_offset_3 o_q, Periodic_3_offset_3 o_r, Periodic_3_offset_3 o_s)`,

which determines the bounded side of the circle defined by
`(p,o_p), (q,o_q)`, and `(r,o_r)` on which the point-offset pair
`(s,o_s)` lies.
\pre `p,q,r`, and `s` are coplanar and `p,q`, and `r` are not collinear, `(p,o_p),(q,o_q),(r,o_r)`, and `(s,o_s)` are coplanar and `(p,o_p),(q,o_q)`, and `(r,o_r)` are not collinear, respectively, and `p`, `q`, `r`, `s`, `t` lie inside the domain.
*/
typedef unspecified_type Coplanar_side_of_bounded_circle_3;

/// @}

/// \name
/// When `is_Gabriel` is used, the traits class must
/// in addition provide the following predicate object:
/// @{

/*!
A predicate object that must provide the function operator

`Bounded_side operator()(Point_3 p, Point_3 q, Point_3 t, Periodic_3_offset_3 o_p, Periodic_3_offset_3 o_q, Periodic_3_offset_3 o_t)`,

which returns the position of the point-offset pair `(t,o_t)`
relative to the sphere that has `(p,o_p)(q,o_q)` as its diameter,

`Bounded_side operator()(Point_3 p, Point_3 q, Point_3 r, Point_3 t, Periodic_3_offset_3 o_p, Periodic_3_offset_3 o_q, Periodic_3_offset_3 o_r, Periodic_3_offset_3 o_q)`,

which returns the position of the point-offset pair `(t,o_t)`
relative to the sphere passing through `(p,o_p), (q,o_q)`, and
`(r,o_r)` and whose center is in the plane defined by these three
point-offset pairs.
*/
typedef unspecified_type Side_of_bounded_sphere_3;

/// @}

/// \name
/// When the dual operations are used, the traits
/// class must in addition provide the following constructor object:
/// @{

/*!
A constructor object that must provide the function operator

`Point_3 operator()(Point_3 p, Point_3 q, Point_3 r, Point_3 s, Periodic_3_offset_3 o_p, Periodic_3_offset_3 o_q, Periodic_3_offset_3 o_r, Periodic_3_offset_3 o_s)`,

which constructs the circumcenter of four point-offset pairs.
\pre `p`, `q`, `r` and `s` as well as `(p,o_p)`, `(q,o_q)`, `(r,o_r)` and `(s,o_s)` must be non coplanar. `p`, `q`, `r`, `s` lie inside the domain.
*/
typedef unspecified_type Construct_circumcenter_3;

/// @}

/// \name Operations
/// The following functions give access to the predicate and construction objects:
/// @{

/*!

*/
Side_of_oriented_sphere_3 side_of_oriented_sphere_3_object();

/*!

*/
Compare_distance_3 compare_distance_3_object();

/// @}

/// \name
/// The following functions must be provided if vertex removal is
/// used; otherwise dummy functions can be provided.
/// @{

/*!

*/
Coplanar_orientation_3 coplanar_3_orientation_3_object();

/*!

*/
Coplanar_side_of_bounded_circle_3
coplanar_side_of_bounded_circle_3_object();

/// @}

/// \name
/// The following function must be provided if the `is_Gabriel`
/// methods of `Periodic_3_Delaunay_triangulation_3` are used;
/// otherwise a dummy function can be provided.
/// @{

/*!

*/
Side_of_bounded_sphere_3 side_of_bounded_sphere_3_object();

/// @}

/// \name
/// The following function must be provided if the methods of
/// `Periodic_3_Delaunay_triangulation_3` returning elements of the
/// Voronoi diagram are used; otherwise a dummy function can be
/// provided:
/// @{

/*!

*/
Construct_circumcenter_3 construct_circumcenter_3_object();

/// @}

}; /* end Periodic_3DelaunayTriangulationTraits_3 */

