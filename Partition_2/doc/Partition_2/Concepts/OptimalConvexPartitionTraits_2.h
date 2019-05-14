/*!
\ingroup PkgPolygonPartitioning2Concepts
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
A segment type 
*/ 
typedef unspecified_type Segment_2; 

/*!
A ray type 
*/ 
typedef unspecified_type Ray_2; 

/*!
A general object type that can be either a point or a segment 
*/ 
typedef unspecified_type Object_2; 

/*!
Function object type that provides 
`Segment_2 operator()(Point_2 p, Point_2 q)`, which constructs and 
returns the segment defined by the points \f$ p\f$ and \f$ q\f$.
*/ 
typedef unspecified_type Construct_segment_2; 

/*!
Function object type that provides 
`Ray_2 operator()(Point_2 p, Point_2 q)`, which constructs and returns 
the ray from point \f$ p\f$ through point \f$ q\f$. 
*/ 
typedef unspecified_type Construct_ray_2; 

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

/*!
Function object type that provides 
`Object_2 operator()(Segment_2 s1, Segment_2 s2)` that returns 
the intersection of two segments (which may be either a segment or 
a point). 
*/ 
typedef unspecified_type Intersect_2; 

/*!
Function object type that provides 
`bool operator()(Segment_2 s1, Object_2 o)` that returns 
`true` if `o` is a segment and assigns the value of `o` 
to `s1`; returns `false` otherwise. 
*/ 
typedef unspecified_type Assign_2; 

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
Construct_segment_2 construct_segment_2_object() const; 

/*!

*/ 
Construct_ray_2 construct_ray_2_object() const; 

/*!

*/ 
Are_strictly_ordered_along_line_2 
are_strictly_ordered_along_line_2_object() const; 

/// @}

}; /* end OptimalConvexPartitionTraits_2 */
