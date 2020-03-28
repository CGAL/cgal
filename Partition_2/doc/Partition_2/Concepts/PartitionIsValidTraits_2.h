/*!
\ingroup PkgPartition2Concepts
\cgalConcept

Requirements of a traits class that is used
by `CGAL::partition_is_valid_2()`, `CGAL::convex_partition_is_valid_2()`,
and `CGAL::y_monotone_partition_is_valid_2()` for testing if a given set of
polygons are nonoverlapping and if their union is a polygon that is the
same as a polygon represented by a given sequence of points. Note that the
traits class for `CGAL::partition_is_valid_2()` may have to satisfy additional
requirements if each partition polygon is to be tested for having a
particular property; see, for example, the descriptions of the
function `CGAL::is_convex_2()`
and the concept `YMonotonePartitionTraits_2` for the additional requirements
for testing for convexity and \f$ y\f$-monotonicity, respectively.

\cgalHasModel `CGAL::Partition_is_valid_traits_2<Traits, PolygonIsValid>`

\sa `CGAL::approx_convex_partition_2()`
\sa `CGAL::greene_approx_convex_partition_2()`
\sa `CGAL::optimal_convex_partition_2()`
\sa `CGAL::y_monotone_partition_2()`

*/

class PartitionIsValidTraits_2 {
public:

/// \name Types
/// @{

/*!
The point type on which the partitioning algorithm operates.
*/
typedef unspecified_type Point_2;

/*!
The polygon type created by the partitioning
function. This type should provide a nested type `Vertex_const_iterator`
that is the type of the non-mutable iterator over the polygon vertices.
*/
typedef unspecified_type Polygon_2;

/*!
A model of the concept `PolygonIsValid`
*/
typedef unspecified_type Is_valid;

/*!

Predicate object type that compares `Point_2`s lexicographically.
Must provide `bool operator()(Point_2 p, Point_2 q)` where `true`
is returned iff \f$ p <_{xy} q\f$.
We have \f$ p<_{xy}q\f$, iff \f$ p_x < q_x\f$ or \f$ p_x = q_x\f$ and \f$ p_y < q_y\f$,
where \f$ p_x\f$ and \f$ p_y\f$ denote the \f$ x\f$ and \f$ y\f$ coordinates of point \f$ p\f$,
respectively.

*/
typedef unspecified_type Less_xy_2;

/*!

Predicate object type that provides
`bool operator()(Point_2 p,Point_2 q,Point_2 r)`, which
returns `true` iff `r` lies to the left of the
oriented line through `p` and `q`.
*/
typedef unspecified_type Left_turn_2;

/*!
Predicate object type that provides
`CGAL::Orientation operator()(Point_2 p, Point_2 q, Point_2 r)` that
returns `CGAL::LEFT_TURN`, if \f$ r\f$ lies to the left of the oriented
line \f$ l\f$ defined by \f$ p\f$ and \f$ q\f$, returns `CGAL::RIGHT_TURN` if \f$ r\f$
lies to the right of \f$ l\f$, and returns `CGAL::COLLINEAR` if \f$ r\f$ lies
on \f$ l\f$.
*/
typedef unspecified_type Orientation_2;

/// @}

/// \name Creation
/// Only a copy constructor is required.
/// @{

/*!

*/
PartitionIsValidTraits_2(PartitionIsValidTraits_2 tr);

/// @}

/// \name Operations
/// The following functions that create instances of the above predicate object types must exist.
/// @{

/*!

*/
Orientation_2 is_valid_object();

/*!

*/
Less_xy_2 less_xy_2_object();

/*!

*/
Left_turn_2 left_turn_2_object();

/*!

*/
Orientation_2 orientation_2_object();

/// @}

}; /* end PartitionIsValidTraits_2 */
