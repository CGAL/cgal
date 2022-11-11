/*!
\ingroup PkgPartition2Concepts
\cgalConcept

Requirements of a traits class to be
used with the function `optimal_convex_partition_2()` that computes
an optimal convex partition of a polygon.

\cgalRefines `PartitionTraits_2`

\cgalHasModel `CGAL::Partition_traits_2<R>`

\sa `CGAL::convex_partition_is_valid_2()`
\sa `CGAL::Partition_is_valid_traits_2<Traits, PolygonIsValid>`

*/

class OptimalConvexPartitionTraits_2 {
public:

/// \name Types
/// In addition to the types listed with the concept
/// `PartitionTraits_2`, the following types are required:
/// @{


/*!
Predicate object type that
determines orderings of \link PartitionTraits_2::Point_2 `Point_2` \endlink 's on a line. Must provide
`bool operator()(Point_2 p, Point_2 q, Point_2 r)` that
returns `true`, iff `q` lies between `p`
and `r` and `p`, `q`, and `r` satisfy the precondition
that they are collinear.
*/
typedef unspecified_type Collinear_are_ordered_along_line_2;

/*!
Predicate object type that
determines orderings of \link PartitionTraits_2::Point_2 `Point_2` \endlink 's. Must provide
`bool operator()(Point_2 p, Point_2 q, Point_2 r)` that
returns `true`, iff the three points are collinear and
`q` lies strictly between `p`
and `r`. Note that `false` should be returned if
`q==p` or `q==r`.
*/
typedef unspecified_type Are_stritcly_ordered_along_line_2;


/// @}

/// \name Creation
/// Only a copy constructor is required.
/// @{

/*!

*/
OptimalConvexPartitionTraits_2(OptimalConvexPartitionTraits_2 tr);

/// @}

/// \name Operations
/// In addition to the functions required by `PartitionTraits_2`, the
/// following functions that create instances of the above function
/// object types must exist:
/// @{

/*!

*/
Collinear_are_ordered_along_line_2 collinear_are_ordered_along_line_2_object() const;


/*!

*/
Are_strictly_ordered_along_line_2
are_strictly_ordered_along_line_2_object() const;

/// @}

}; /* end OptimalConvexPartitionTraits_2 */
