/*!
\ingroup PkgPartition2Concepts
\cgalConcept

The polygon partitioning functions are each parameterized by a traits class
that defines the primitives used in the algorithms. Many requirements are
common
to all traits classes. The concept `PartitionTraits_2` defines this common set of
requirements.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Partition_traits_2<R>}
\cgalHasModelsEnd

\sa `CGAL::approx_convex_partition_2()`
\sa `CGAL::greene_approx_convex_partition_2()`
\sa `CGAL::optimal_convex_partition_2()`
\sa `CGAL::y_monotone_partition_2()`

*/

class PartitionTraits_2 {
public:

/// \name Types
/// @{

/*!
The point type on which the partitioning algorithm operates.
*/
typedef unspecified_type Point_2;

/*!
The polygon type to be created by the partitioning
algorithm. For testing the validity postcondition of the partition, this
type should provide a nested type `Vertex_const_iterator` that is the
type of the iterator over the polygon vertices and member functions
`Vertex_const_iterator vertices_begin()` and
`Vertex_const_iterator vertices_end()`.
*/
typedef unspecified_type Polygon_2;

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

Same as `Less_xy_2` with the roles of \f$ x\f$ and \f$ y\f$ interchanged.
*/
typedef unspecified_type Less_yx_2;

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

/*!
Predicate object type that provides
`CGAL::Comparison_result operator()(Point_2 p, Point_2 q)` to compare
the \f$ y\f$ values of two points. The operator must return
`CGAL::SMALLER` if \f$ p_y < q_y\f$, `CGAL::LARGER` if \f$ p_y > q_y\f$ and
`CGAL::EQUAL` if \f$ p_y = q_y\f$.
*/
typedef unspecified_type Compare_y_2;

/*!
The same as `Compare_y_2`, except that \f$ x\f$
coordinates are compared instead of \f$ y\f$.
*/
typedef unspecified_type Compare_x_2;

/// @}

/// \name Creation
/// A copy constructor and default constructor are required.
/// @{

/*!

*/
PartitionTraits_2();

/*!

*/
PartitionTraits_2(PartitionTraits_2 tr);

/// @}

/// \name Operations
/// The following functions that create instances of the above predicate object types must exist.
/// @{

/*!

*/
Less_yx_2 less_yx_2_object() const;

/*!

*/
Less_xy_2 less_xy_2_object() const;

/*!

*/
Left_turn_2 left_turn_2_object() const;

/*!

*/
Orientation_2 orientation_2_object() const;

/*!

*/
Compare_y_2 compare_y_2_object() const;

/*!

*/
Compare_x_2 compare_x_2_object() const;

/// @}

}; /* end PartitionTraits_2 */
