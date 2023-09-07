
/*!
\ingroup PkgTriangulation3Concepts
\cgalConcept

The concept `RegularTriangulationTraits_3` is the first template parameter of the class
`CGAL::Regular_triangulation_3`. It defines the geometric objects (points,
segments...) forming the triangulation together with a few geometric
predicates and constructions on these objects.

We use here the same notation as in Section \ref
Triangulation3secclassRegulartriangulation. To simplify notation, \f$
p\f$ will often denote in the sequel either the point \f$ p\in\mathbb{R}^3\f$
or the weighted point \f$ {p}^{(w)}=(p,w_p)\f$.

\cgalRefines{TriangulationTraits_3}

\cgalHasModelsBegin
\cgalHasModelsBare{All models of the \cgal concept `Kernel`}
\cgalHasModelsEnd

\sa `CGAL::Regular_triangulation_3`

In addition to the requirements described for the traits class of
  `CGAL::Triangulation_3`, the geometric traits class of
  `CGAL::Regular_triangulation_3` must fulfill the following requirements.

*/
class RegularTriangulationTraits_3 {
public:

/// \name Types
/// @{

/*!
The line type.
*/
typedef unspecified_type Line_3;

/*!
The object type.
*/
typedef unspecified_type Object_3;

/*!
The plane type.
*/
typedef unspecified_type Plane_3;

/*!
The ray type.
*/
typedef unspecified_type Ray_3;

/*!
The weighted point type. It has to be a model of the concept `Kernel::WeightedPoint_3`.

\note The unweighted point type `Point_3` is requested by the concept
`TriangulationTraits_3`, which this concept refines.
*/
typedef unspecified_type Weighted_point_3;

/*!
A predicate object,
model of `Kernel::PowerSideOfOrientedPowerSphere_3`,
that must provide the following function operators:

`Oriented_side operator()( Weighted_point_3 p,                          Weighted_point_3 q,                          Weighted_point_3 r,                          Weighted_point_3 s,                          Weighted_point_3 t)`,

which performs the following:

Let \f$ {z(p,q,r,s)}^{(w)}\f$ be the power sphere of the weighted points
\f$ (p,q,r,s)\f$. Returns

- `ON_ORIENTED_BOUNDARY` if `t` is orthogonal to
  \f$ {z(p,q,r,s)}^{(w)}\f$,

- `ON_NEGATIVE_SIDE` if `t` lies outside the oriented sphere of
  center \f$ z(p,q,r,s)\f$ and radius \f$ \sqrt{ w_{z(p,q,r,s)}^2 + w_t^2 }\f$
  (which is equivalent to \f$ \Pi({t}^{(w)},{z(p,q,r,s)}^{(w)}) >0\f$),

- `ON_POSITIVE_SIDE` if `t` lies inside this oriented sphere.

\pre `p, q, r, s` are not coplanar.
Note that with this definition, if all the points have a weight equal
to 0, then
`power_side_of_oriented_power_sphere_3(p,q,r,s,t)` = `side_of_oriented_sphere(p,q,r,s,t)`.

<HR WIDTH=50%>

`Oriented_side operator()( Weighted_point_3 p,                          Weighted_point_3 q,                          Weighted_point_3 r,                          Weighted_point_3 t)`,

which has a
definition analogous to the previous method, for coplanar points,
with the power circle \f$ {z(p,q,r)}^{(w)}\f$.
\pre `p, q, r` are not collinear and `p, q, r, t` are coplanar.
If all the points have a weight equal to 0, then
`power_side_of_oriented_power_sphere_3(p,q,r,t)` = `side_of_oriented_circle(p,q,r,t)`.

<HR WIDTH=50%>

`Oriented_side operator()( Weighted_point_3 p,                          Weighted_point_3 q,                          Weighted_point_3 t)`,

which is the same for collinear points, where \f$ {z(p,q)}^{(w)}\f$ is the
power segment of `p` and `q`.
\pre `p` and `q` have different bare points, and `p, q, t` are collinear.
If all points have a weight equal to 0, then
`power_side_of_oriented_power_sphere_3(p,q,t)` gives the same answer as the kernel predicate
`s(p,q).has_on(t)` would give, where `s(p,q)` denotes the
segment with endpoints `p` and `q`.

<HR WIDTH=50%>

`Oriented_side operator()( Weighted_point_3 p, Weighted_point_3 q)`,

which is the same for equal bare points, then it returns the comparison of the weights
(`ON_POSITIVE_SIDE` when `q` is heavier than `p`).
\pre `p` and `q` have equal bare points.

*/
typedef unspecified_type Power_side_of_oriented_power_sphere_3;


/*!
A predicate object,
model of `Kernel::ComparePowerDistance_3`,
that must provide the function operator

`Comparison_result operator()(Point_3 p, Weighted_point_3 q, Weighted_point_3 r)`,

which compares the power distance between `p` and `q`
to the power distance
between `p` and `r`.

\note This predicate is required if a call to
`nearest_power_vertex()` or `nearest_power_vertex_in_cell()` is
issued.
*/
typedef unspecified_type Compare_power_distance_3;

/*!
A constructor type,
model of `Kernel::ConstructPoint_3`.
The `operator()` extracts the bare point from a weighted point.

`Point_3 operator() ( Weighted_point_3 p);`
*/
typedef unspecified_type Construct_point_3;

/*!
A constructor type,
model of `Kernel::ConstructWeightedCircumcenter_3`.
The `operator()` constructs the bare point
which is the center of the smallest orthogonal sphere to the input
weighted points.

`Point_3 operator() ( Weighted_point_3 p, Weighted_point_3 q, Weighted_point_3 r, Weighted_point_3 s);`


\note Only required when the dual operations are used.
*/
typedef unspecified_type Construct_weighted_circumcenter_3;

/*!
A constructor object that must provide the function operators

`Object_3 operator()(Point_3 p)`,

`Object_3 operator()(Segment_3 s)` and

`Object_3 operator()(Ray_3 r)`

that construct an object respectively from a point, a segment and a ray.

\note Only required when the dual operations are used.
*/
typedef unspecified_type Construct_object_3;

/*!
A constructor object that must provide the function operator

`Line_3 operator()(Plane_3 pl, Point_3 p)`,

which constructs the line perpendicular to `pl` passing through `p`.

\note Only required when the dual operations are used.
*/
typedef unspecified_type Construct_perpendicular_line_3;

/*!
A constructor object that must provide the function operator

`Plane_3 operator()(Point_3 p, Point_3 q, Point_3 r)`,

which constructs the plane passing through `p`, `q` and `r`.
\pre `p`, `q` and `r` are non collinear.

\note Only required when the dual operations are used.
*/
typedef unspecified_type Construct_plane_3;

/*!
A constructor object that must provide the function operator

`Ray_3 operator()(Point_3 p, Line_3 l)`,

which constructs the ray starting at `p` with direction given by `l`.

\note Only required when the dual operations are used.
*/
typedef unspecified_type Construct_ray_3;

/// @}

/// \name
/// When `is_Gabriel` functions are used, the traits class must
/// in addition provide the following predicate object:
/// @{

/*!
A predicate object that must provide the function operators

`Bounded_side operator()(Weighted_point_3 p, Weighted_point_3 t)`,

which returns the sign of the power test of `t` with respect to the smallest
sphere orthogonal to `p` (which is the sphere with center `p` and squared
radius `-w_p` with `w_p` the weight of `p`),

`Bounded_side operator()(Weighted_point_3 p, Weighted_point_3 q, Weighted_point_3 t)`,

which returns the sign of the power test of `t` with respect to the smallest
sphere orthogonal to `p` and `q`,

`Bounded_side operator()(Weighted_point_3 p, Weighted_point_3 q, Weighted_point_3 r, Weighted_point_3 t)`,

which returns the sign of the power test of `t` with respect to the smallest
sphere orthogonal to `p`, `q`, and `r`.
*/
typedef unspecified_type Power_side_of_bounded_power_sphere_3;

/// @}

/// \name Operations
/// @{

/*!

*/
Power_side_of_oriented_power_sphere_3 power_side_of_oriented_power_sphere_3_object();

/*!

*/
Compare_power_distance_3 compare_power_distance_3_object();

/*!

*/
Construct_point_3 construct_point_3_object();

/// @}

/*! \name
The following functions must be provided only if the member functions of
`CGAL::Regular_triangulation_3` returning elements of the dual diagram are called:
*/
/// @{


Construct_weighted_circumcenter_3 construct_weighted_circumcenter_3_object();

/*!

*/
Construct_object_3 construct_object_3_object();

/*!

*/
Construct_perpendicular_line_3 construct_perpendicular_line_object();

/*!

*/
Construct_plane_3 construct_plane_3_object();

/*!

*/
Construct_ray_3 construct_ray_3_object();

/// @}

}; /* end RegularTriangulationTraits_3 */

