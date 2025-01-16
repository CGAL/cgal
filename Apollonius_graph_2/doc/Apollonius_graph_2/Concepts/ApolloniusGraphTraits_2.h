
/*!
\ingroup PkgApolloniusGraph2Concepts
\cgalConcept

\cgalRefines{TriangulationTraits_2}

The concept `ApolloniusGraphTraits_2` provides the traits
requirements for the `Apollonius_graph_2` class. In particular,
it provides a type `Site_2`, which must be a model of the concept
`ApolloniusSite_2`. It also provides
constructions for sites and several function object
types for the predicates.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Apollonius_graph_traits_2<K,Method_tag>}
\cgalHasModels{CGAL::Apollonius_graph_filtered_traits_2<CK,CM,EK,EM,FK,FM>}
\cgalHasModelsEnd

\sa `CGAL::Apollonius_graph_2<Gt,Agds>`
\sa `CGAL::Apollonius_graph_traits_2<K,Method_tag>`
\sa `CGAL::Apollonius_graph_filtered_traits_2<CK,CM,EK,EM,FK,FM>`

*/

class ApolloniusGraphTraits_2 {
public:

/// \name Types
/// @{

/*!
A type for a point.
*/
typedef unspecified_type Point_2;

/*!
A type for an Apollonius site. Must be a model
of the concept `ApolloniusSite_2`.
*/
typedef unspecified_type Site_2;

/*!
A type for a line. Only required if access to
the dual of the Apollonius graph is required or if the primal
or dual diagram are inserted in a stream.
*/
typedef unspecified_type Line_2;

/*!
A type for a ray. Only required if access to
the dual of the Apollonius graph is required or if the primal
or dual diagram are inserted in a stream.
*/
typedef unspecified_type Ray_2;

/*!
A type for a segment. Only required if access to
the dual of the Apollonius graph is required or if the primal
or dual diagram are inserted in a stream.
*/
typedef unspecified_type Segment_2;

/*!
A type representing different types of objects
in two dimensions, namely: `Point_2`, `Site_2`,
`Line_2`, `Ray_2` and `Segment_2`.
*/
typedef unspecified_type Object_2;

/*!
A type for the field number type of sites.
*/
typedef unspecified_type FT;

/*!
A type for the ring number type of sites.
*/
typedef unspecified_type RT;

/*!
Must provide `template <class T> bool operator() ( T& t,
Object_2 o)` which assigns `o` to `t` if `o` was
constructed from an object of type `T`. Returns
`true`, if the assignment was possible.
*/
typedef unspecified_type Assign_2;

/*!
Must provide `template <class T>
Object_2 operator()( T t)` that constructs an object of type
`Object_2` that contains `t` and returns it.
*/
typedef unspecified_type Construct_object_2;

/*!

A constructor for a point of the Apollonius diagram equidistant
from three sites. Must provide
`Point_2 operator()(Site_2 s1, Site_2 s2, Site_2 s3)`, which
constructs a point equidistant from the sites `s1`, `s2` and
`s3`.

*/
typedef unspecified_type Construct_Apollonius_vertex_2;

/*!
A constructor for
a dual Apollonius site (a site whose center is a
vertex of the Apollonius diagram and its weight is the common
distance of its center from the three defining sites).
Must provide `Site_2 operator()(Site_2 s1,
Site_2 s2, Site_2 s3)`, which constructs a
dual site whose center \f$ c\f$ is equidistant from `s1`, `s2` and
`s3`, and its weight is equal to the (signed) distance of \f$ c\f$
from `s1` (or `s2` or `s3`).

Must also provide `Line_2 operator()(Site_2 s1, Site_2 s2)`, which
constructs a line bitangent to `s1` and `s2`. This line is the
dual site of `s1`, `s2` and the site at infinity; it can be
viewed as a dual Apollonius site whose center is at infinity
and its weight is infinite.

*/
typedef unspecified_type Construct_Apollonius_site_2;

/*!
A predicate object type. Must
provide `Comparison_result operator()(Site_2 s1,
Site_2 s2)`, which compares the \f$ x\f$-coordinates of the centers of
`s1` and `s2`.
*/
typedef unspecified_type Compare_x_2;

/*!
A predicate object type. Must
provide `Comparison_result operator()(Site_2 s1,
Site_2 s2)`, which compares the \f$ y\f$-coordinates of the centers of
`s1` and `s2`.
*/
typedef unspecified_type Compare_y_2;

/*!
A predicate object type. Must
provide `Comparison_result operator()(Site_2 s1,
Site_2 s2)`, which compares the weights of `s1`
and `s2`.
*/
typedef unspecified_type Compare_weight_2;

/*!
A predicate object type. Must
provide `Orientation operator()(Site_2 s1,
Site_2 s2, Site_2 s3)`, which performs the
usual orientation test for the centers of the three sites
`s1`, `s2` and `s3`.

Must also provide `Orientation operator()(Site_2 s1, Site_2 s2,
Site_2 s3, Site_2 p1, Site_2 p2)`,
which performs the usual orientation test for the Apollonius vertex of
`s1`, `s2`, `s3` and the centers of `p1` and
`p2`.
\pre the Apollonius vertex of `s1`, `s2` and `s3` must exist.
*/
typedef unspecified_type Orientation_2;

/*!
A predicate object type. Must
provide `bool operator()(Site_2 s1,
Site_2 s2)`, which returns `true` if the circle
corresponding to `s2` is contained in the closure of the disk
corresponding to `s1`, `false` otherwise.
*/
typedef unspecified_type Is_hidden_2;

/*!
A predicate object type.
Must provide `Oriented_side operator()(Site_2 s1,
Site_2 s2, Point_2 p)`, which returns
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
Must provide `Sign operator()(Site_2 s1, Site_2
s2, Site_2 s3, Site_2 q)`, which
returns the sign of the distance of `q` from the dual Apollonius
site of `s1`, `s2`, `s3`.
\pre the dual Apollonius site of `s1`, `s2`, `s3` must exist.

Must also provide `Sign operator()(Site_2 s1,
Site_2 s2, Site_2 q)`, which returns the sign of the distance of
`q` from the bitangent line of `s1`, `s2` (a degenerate
dual Apollonius site, with its center at infinity).
*/
typedef unspecified_type Vertex_conflict_2;

/*!
A predicate object type.
Must provide `bool operator()(Site_2 s1, Site_2 s2, Site_2 s3, Site_2 s4, Site_2 q, bool b)`.
The sites `s1`, `s2`, `s3` and `s4` define an Apollonius edge that lies on the
bisector of `s1` and `s2` and has as endpoints the Apollonius
vertices defined by the triplets `s1`, `s2`, `s3` and
`s1`, `s4` and `s2`. The Boolean `b` denotes if the
two Apollonius vertices are in conflict with the site
`q` (in which case `b` should be `true`, otherwise
`false`). In case that `b` is `true`, the predicate
returns `true` if and only if the entire Apollonius edge is in
conflict with `q`. If `b` is `false`, the predicate returns
`false` if and only if `q` is not in conflict with the
Apollonius edge.
\pre the Apollonius vertices of `s1`, `s2`, `s3`, and `s1`, `s4`, `s2` must exist.

Must also provide `bool operator()(Site_2 s1,
Site_2 s2, Site_2 s3, Site_2 q, bool b)`. The
sites `s1`, `s2`, `s3` and the site at infinity
\f$ s_\infty\f$ define an Apollonius edge that lies on the bisector of
`s1` and `s2` and has as endpoints the Apollonius vertices
defined by the triplets `s1`, `s2`, `s3` and `s1`,
\f$ s_\infty\f$ and `s2` (the second Apollonius vertex is actually at
infinity). The Boolean `b` denotes if the two Apollonius vertices
are in conflict with the site `q` (in which case `b`
should be `true`, otherwise `false`).
In case that `b` is `true`, the predicate
returns `true` if and only if the entire Apollonius edge is in
conflict with `q`. If `b` is `false`, the predicate returns
`false` if and only if `q` is not in conflict with the
Apollonius edge.
\pre the Apollonius vertex of `s1`, `s2`, `s3` must exist.

Must finally provide `bool operator()(Site_2 s1,
Site_2 s2, Site_2 q, bool b)`. The
sites `s1`, `s2` and the site at infinity
\f$ s_\infty\f$ define an Apollonius edge that lies on the bisector of
`s1` and `s2` and has as endpoints the Apollonius vertices
defined by the triplets `s1`, `s2`, \f$ s_\infty\f$ and `s1`,
\f$ s_\infty\f$ and `s2` (both Apollonius vertices are actually at
infinity). The Boolean `b` denotes if the two Apollonius vertices
are in conflict with the site `q` (in which case `b`
should be `true`, otherwise `false`).
In case that `b` is `true`, the predicate
returns `true` if and only if the entire Apollonius edge is in
conflict with `q`. If `b` is `false`, the predicate returns
`false` if and only if `q` is not in conflict with the
Apollonius edge.
*/
typedef unspecified_type Finite_edge_interior_conflict_2;

/*!
A predicate
object type. Must provide `bool operator()(Site_2 s1,
Site_2 s2, Site_2 s3, Site_2 q, bool b)`. The
sites \f$ s_\infty\f$, `s1`, `s2` and `s3` define an
Apollonius edge that lies on the bisector of \f$ s_\infty\f$ and `s1`
and has as endpoints the Apollonius vertices defined by the triplets
\f$ s_\infty\f$, `s1`, `s2` and \f$ s_\infty\f$, `s3` and
`s1`. The Boolean `b` denotes if the two Apollonius vertices
are in conflict with the site `q` (in which case `b`
should be `true`, otherwise `false`.
In case that `b` is `true`, the predicate
returns `true` if and only if the entire Apollonius edge is in
conflict with `q`. If `b` is `false`, the predicate returns
`false` if and only if `q` is not in conflict with the
Apollonius edge.
*/
typedef unspecified_type Infinite_edge_interior_conflict_2;

/*!
A predicate object type.
Must provide `bool operator()(Site_2 s1, Site_2
s2, Site_2 s3, Site_2 s4)`. It returns `true` if
the Apollonius edge defined by `s1`, `s2`, `s3` and
`s4` is degenerate, `false` otherwise. An Apollonius edge is
called degenerate if its two endpoints coincide.
\pre the Apollonius vertices of `s1`, `s2`, `s3`, and `s1`, `s4`, `s2` must exist.
*/
typedef unspecified_type Is_degenerate_edge_2;

/// @}

/// \name Creation
/// @{

/*!
Default constructor.
*/
ApolloniusGraphTraits_2();

/*!
Copy constructor.
*/
ApolloniusGraphTraits_2(ApolloniusGraphTraits_2 other);

/*!
Assignment operator.
*/
ApolloniusGraphTraits_2 operator=(ApolloniusGraphTraits_2 other);

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
Compare_weight_2 compare_weight_2_object();

/*!

*/
Orientation_2 orientation_2_object();

/*!

*/
Is_hidden_2 is_hidden_2_object();

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
Is_degenerate_edge_2 is_degenerate_edge_2_object();

/// @}

/// \name Access to constructor objects
/// @{

/*!

*/
Construct_object_2
construct_object_2_object();

/*!

*/
Construct_Apollonius_vertex_2
construct_Apollonius_vertex_2_object();

/*!

*/
Construct_Apollonius_site_2
construct_Apollonius_site_2_object();

/// @}

/// \name Access to other objects
/// @{

/*!

*/
Assign_2 assign_2_object();

/// @}

}; /* end ApolloniusGraphTraits_2 */

