
/*!
\ingroup PkgSegmentDelaunayGraph2Concepts

\cgalConcept

\cgalRefines `TriangulationTraits_2`

The concept `SegmentDelaunayGraphTraits_2` provides the traits
requirements for the `CGAL::Segment_Delaunay_graph_2<Gt,St,DS>` and
`CGAL::Segment_Delaunay_graph_hierarchy_2<Gt,St,STag,DS>` classes.
In particular, it provides a type `Site_2`, which must be a model of
the concept `SegmentDelaunayGraphSite_2`. It also provides
constructions for sites and several function object
types for the predicates.

\cgalHasModel `CGAL::Segment_Delaunay_graph_traits_2<K,MTag>`
\cgalHasModel `CGAL::Segment_Delaunay_graph_traits_without_intersections_2<K,MTag>`
\cgalHasModel `CGAL::Segment_Delaunay_graph_filtered_traits_2<CK,CM,EK,EM,FK,FM>`
\cgalHasModel `CGAL::Segment_Delaunay_graph_filtered_traits_without_intersections_2<CK,CM,EK,EM,FK,FM>`

\sa `SegmentDelaunayGraphSite_2`
\sa `CGAL::Segment_Delaunay_graph_2<Gt,St,DS>`
\sa `CGAL::Segment_Delaunay_graph_hierarchy_2<Gt,St,STag,DS>`
\sa `CGAL::Segment_Delaunay_graph_traits_2<K,MTag>`
\sa `CGAL::Segment_Delaunay_graph_traits_without_intersections_2<K,MTag>`
\sa `CGAL::Segment_Delaunay_graph_filtered_traits_2<CK,CM,EK,EM,FK,FM>`
\sa `CGAL::Segment_Delaunay_graph_filtered_traits_without_intersections_2<CK,CM,EK,EM,FK,FM>`

*/

class SegmentDelaunayGraphTraits_2 {
public:

/// \name Types
/// @{

/*!
Indicates or not whether the
intersecting segments are to be supported. The tag must either be
`CGAL::Tag_true` or `CGAL::Tag_false`.
*/
typedef unspecified_type Intersections_tag;

/*!
A type for a site of the segment Delaunay
graph. Must be a model of the concept
`SegmentDelaunayGraphSite_2`.
*/
typedef unspecified_type Site_2;

/*!
A type for a point.
*/
typedef unspecified_type Point_2;

/*!
A type for a line. Only required if the segment
Delaunay graph is inserted in a stream.
*/
typedef unspecified_type Line_2;

/*!
A type for a ray. Only required if the segment
Delaunay graph is inserted in a stream.
*/
typedef unspecified_type Ray_2;

/*!
A type for a segment. Only required if
if the segment Delaunay graph is inserted in a stream.
*/
typedef unspecified_type Segment_2;

/*!
A type for the field number type of sites, points, etc...
*/
typedef unspecified_type FT;

/*!
A type for the ring number type of sites, points, etc.
*/
typedef unspecified_type RT;

/*!
An enumeration type that indicates the
type of the arrangement of two sites. The possible values are
`DISJOINT`, `IDENTICAL`, `CROSSING`,
`TOUCHING_1`, `TOUCHING_2`, `TOUCHING_11`,
`TOUCHING_12`, `TOUCHING_21`, `TOUCHING_22`,
`OVERLAPPING_11`, `OVERLAPPING_12`, `OVERLAPPING_21`,
`OVERLAPPING_22`, `INTERIOR`, `INTERIOR_1`,
`INTERIOR_2`, `TOUCHING_11_INTERIOR_1`,
`TOUCHING_11_INTERIOR_2`, `TOUCHING_12_INTERIOR_1`,
`TOUCHING_12_INTERIOR_2`, `TOUCHING_21_INTERIOR_1`,
`TOUCHING_21_INTERIOR_2`, `TOUCHING_22_INTERIOR_1`,
`TOUCHING_22_INTERIOR_2`. A detailed description of the meaning
of these values is shown the end of the reference manual for this
concept.
*/
typedef unspecified_type Arrangement_type;

/*!
A type representing different types of objects
in two dimensions, namely: `Point_2`, `Site_2`,
`Line_2`, `Ray_2` and `Segment_2`.
*/
typedef unspecified_type Object_2;

/*!
Must provide `template <class T> bool operator() ( T& t, Object_2 o)` which assigns `o` to `t` if `o` was
constructed from an object of type `T`. Returns
`true`, if the assignment was possible.
*/
typedef unspecified_type Assign_2;

/*!
Must provide `template <class T> Object_2 operator()( T t)` that constructs an object of type
`Object_2` that contains `t` and returns it.
*/
typedef unspecified_type Construct_object_2;

/*!

A constructor for a point of the segment Voronoi diagram equidistant
from three sites. Must provide
`Point_2 operator()(Site_2 s1, Site_2 s2, Site_2 s3)`, which
constructs a point equidistant from the sites `s1`, `s2` and
`s3`.

*/
typedef unspecified_type Construct_svd_vertex_2;

/*!
A predicate object type. Must
provide `Comparison_result operator()(Site_2 s1, Site_2 s2)`, which compares the \f$ x\f$-coordinates of the points
represented by the sites `s1` and `s2`.
\pre `s1` and `s2` must be points.
*/
typedef unspecified_type Compare_x_2;

/*!
A predicate object type. Must
provide `Comparison_result operator()(Site_2 s1, Site_2 s2)`, which compares the \f$ y\f$-coordinates of the points
represented by the sites `s1` and `s2`.
\pre `s1` and `s2` must be points.
*/
typedef unspecified_type Compare_y_2;

/*!
A predicate object type. Must
provide `bool operator()(Point_2 p1, Point_2 p2)`, which returns `true` if `p1.x() < p2.x()`.
*/
typedef unspecified_type Less_x_2;

/*!
A predicate object type. Must
provide `bool operator()(Point_2 p1, Point_2 p2)`, which returns `true` if `p1.y() < p2.y()`.
*/
typedef unspecified_type Less_y_2;

/*!
A predicate object type. Must
provide `Orientation operator()(Site_2 s1, Site_2 s2, Site_2 s3)`, which performs the
usual orientation test for three points.
`s1`, `s2` and `s3`.
\pre the sites `s1`, `s2` and `s3` must be points.
*/
typedef unspecified_type Orientation_2;

/*!
A predicate object type. Must provide
`bool operator()(Site_2 s1, Site_2 s2)`, which determines is the
points represented by the sites `s1` and `s2` are identical.
\pre `s1` and `s2` must be points.
*/
typedef unspecified_type Equal_2;

/*!
A predicate object type. Must provide
`bool operator()(Site_2 s1, Site_2 s2)`, which determines is the
segments represented by the sites `s1` and `s2` are
parallel.
\pre `s1` and `s2` must be segments.
*/
typedef unspecified_type Are_parallel_2;

/*!
A predicate object type.
Must provide `Oriented_side operator()(Site_2 s1, Site_2 s2, Point_2 p)`, which returns
the oriented side of the bisector of `s1` and `s2` that
contains `p`. Returns `ON_POSITIVE_SIDE` if `p` lies in
the half-space of `s1` (i.e., `p` is closer to `s1` than
`s2`); returns `ON_NEGATIVE_SIDE` if `p` lies in the
half-space of `s2`; returns `ON_ORIENTED_BOUNDARY` if `p`
lies on the bisector of `s1` and `s2`.
*/
typedef unspecified_type Oriented_side_of_bisector_2;

/*!
A predicate object type.
Must provide `Sign operator()(Site_2 s1, Site_2 s2, Site_2 s3, Site_2 q)`, which
returns the sign of the distance of `q` from the Voronoi circle
of `s1`, `s2`, `s3` (the Voronoi circle of three sites
`s1`, `s2`, `s3` is a circle co-tangent to all three
sites, that touches them in that order as we walk on its circumference
in the counter-clockwise sense).
\pre the Voronoi circle of `s1`, `s2`, `s3` must exist.

Must also provide `Sign operator()(Site_2 s1, Site_2 s2, Site_2 q)`, which returns the sign of the distance of
`q` from the bitangent line of `s1`, `s2` (a degenerate
Voronoi circle, with its center at infinity).
*/
typedef unspecified_type Vertex_conflict_2;

/*!
A predicate object
type. Must provide `bool operator()(Site_2 s1, Site_2 s2, Site_2 s3, Site_2 s4, Site_2 q, Sign sgn)`. The sites `s1`, `s2`,
`s3` and `s4` define a Voronoi edge that lies on the
bisector of `s1` and `s2` and has as endpoints the Voronoi
vertices defined by the triplets `s1`, `s2`, `s3` and
`s1`, `s4` and `s2`. The sign `sgn` is the common sign
of the distance of the site `q` from the Voronoi circle of the
triplets `s1`, `s2`, `s3` and `s1`, `s4` and
`s2`. In case that `sgn` is equal to `NEGATIVE`, the
predicate returns `true` if and only if the entire Voronoi edge is
in conflict with `q`. If `sgn` is equal to `POSITIVE` or
`ZERO`, the predicate returns `false` if and only if `q`
is not in conflict with the Voronoi edge.
\pre the Voronoi vertices of `s1`, `s2`, `s3`, and `s1`, `s4`, `s2` must exist.

Must also provide `bool operator()(Site_2 s1, Site_2 s2, Site_2 s3, Site_2 q, Sign sgn)`. The
sites `s1`, `s2`, `s3` and the site at infinity
\f$ s_\infty\f$ define a Voronoi edge that lies on the bisector of
`s1` and `s2` and has as endpoints the Voronoi vertices
\f$ v_{123}\f$ and \f$ v_{1\infty{2}}\f$ defined by the triplets `s1`,
`s2`, `s3` and `s1`, \f$ s_\infty\f$ and `s2` (the second
vertex is actually at infinity). The sign `sgn` is the common sign
of the distance of the site `q` from the two Voronoi circles
centered at the Voronoi vertices \f$ v_{123}\f$ and \f$ v_{1\infty{2}}\f$.
In case that `sgn` is `NEGATIVE`, the predicate
returns `true` if and only if the entire Voronoi edge is in
conflict with `q`. If `sgn` is `POSITIVE` or `ZERO`,
the predicate returns `false` if and only if `q` is not in
conflict with the Voronoi edge.
\pre the Voronoi vertex \f$ v_{123}\f$ of `s1`, `s2`, `s3` must exist.

Must finally provide `bool operator()(Site_2 s1, Site_2 s2, Site_2 q, Sign sgn)`. The
sites `s1`, `s2` and the site at infinity
\f$ s_\infty\f$ define a Voronoi edge that lies on the bisector of
\f$ v_{12\infty}\f$ and \f$ v_{1\infty{}2}\f$
`s1` and `s2` and has as endpoints the Voronoi vertices
defined by the triplets `s1`, `s2`, \f$ s_\infty\f$ and `s1`,
\f$ s_\infty\f$ and `s2` (both vertices are actually at
infinity). The sign `sgn` denotes the common sign of the distance
of the site `q` from the Voronoi circles centered at
\f$ v_{12\infty}\f$ and \f$ v_{1\infty{}2}\f$.
If `sgn` is `NEGATIVE`, the predicate
returns `true` if and only if the entire Voronoi edge is in
conflict with `q`. If `POSITIVE` or `ZERO` is `false`,
the predicate returns `false` if and only if `q` is not in
conflict with the Voronoi edge.
*/
typedef unspecified_type Finite_edge_interior_conflict_2;

/*!
A predicate
object type. Must provide `bool operator()(Site_2 s1, Site_2 s2, Site_2 s3, Site_2 q, Sign sgn)`. The
sites \f$ s_\infty\f$, `s1`, `s2` and `s3` define a
Voronoi edge that lies on the bisector of \f$ s_\infty\f$ and `s1`
and has as endpoints the Voronoi vertices \f$ v_{\infty{}12}\f$ and
\f$ v_{\infty{}31}\f$ defined by the triplets
\f$ s_\infty\f$, `s1`, `s2` and \f$ s_\infty\f$, `s3` and
`s1`. The sign `sgn` is the common sign of the distances of
`q` from the Voronoi circles centered at the vertices
\f$ v_{\infty{}12}\f$ and \f$ v_{\infty{}31}\f$. If `sgn` is `NEGATIVE`,
the predicate returns `true` if and only if the entire Voronoi
edge is in conflict with `q`. If `sgn` is `POSITIVE` or
`ZERO`, the predicate returns `false` if and only if `q`
is not in conflict with the Voronoi edge.
*/
typedef unspecified_type Infinite_edge_interior_conflict_2;

/*!
A predicate object type. Must provide
`Oriented_side operator()(Site_1 s1, Site_2 s2, Site_2 s3, Site_2 s, Site_2 p)`. Determines the oriented side of the line \f$ \ell\f$
that contains the point site `p`, where \f$ \ell\f$ is the line
that passes through the Voronoi vertex of the sites `s1`,
`s2`, `s3` and is perpendicular to the segment site `s`.
\pre `s` must be a segment and `p` must be a point.
*/
typedef unspecified_type Oriented_side_2;

/*!
A predicate object type. Must
provide `Arrangement_type operator()(Site_2 s1, Site_2 s2)` that
returns the type of the arrangement of the two sites `s1` and
`s2`.
*/
typedef unspecified_type Arrangement_type_2;

/// @}

/// \name Access to predicate objects
/// @{

/*!

*/
Compare_x_2 compare_x_2_object();

/*!

*/
Compare_y_2 compare_y_2_object();

/*!

*/
Less_x_2 less_x_2_object();

/*!

*/
Less_y_2 less_y_2_object();

/*!

*/
Orientation_2 orientation_2_object();

/*!

*/
Equal_2 equal_2_object();

/*!

*/
Are_parallel_2 are_parallel_2_object();

/*!

*/
Oriented_side_of_bisector_2
oriented_side_of_bisector_test_2_object();

/*!

*/
Vertex_conflict_2 vertex_conflict_2_object();

/*!

*/
Finite_edge_interior_conflict_2
finite_edge_interior_conflict_2_object();

/*!

*/
Infinite_edge_interior_conflict_2
infinite_edge_interior_conflict_2_object();

/*!

*/
Oriented_side_2 oriented_side_2_object();

/*!

*/
Arrangement_type_2 arrangement_type_2_object();

/// @}

/// \name Access to contructor objects
/// @{

/*!

*/
Construct_object_2
construct_object_2_object();

/*!

*/
Construct_svd_vertex_2
construct_svd_vertex_2_object();

/// @}

/// \name Access to other objects
/// @{

/*!

*/
Assign_2 assign_2_object();

/// @}

}; /* end SegmentDelaunayGraphTraits_2 */

