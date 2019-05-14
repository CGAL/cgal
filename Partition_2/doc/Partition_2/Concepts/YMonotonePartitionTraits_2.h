/*!
\ingroup PkgPolygonPartitioning2Concepts
\cgalConcept

Requirements of a traits class to be 
used with the function `y_monotone_partition_2()`. 

\cgalRefines `PartitionTraits_2` 

\cgalHasModel `CGAL::Partition_traits_2<R>` 

*/

class YMonotonePartitionTraits_2 {
public:

/// \name Types 
/// In addition to the types defined for the concept `PartitionTraits_2`, the following types are also required:
/// @{

/*!

*/ 
typedef unspecified_type Line_2; 

/*!
Predicate object type that provides 
`CGAL::Comparision_result operator()(Point_2 p, Line_2 h)` to compare 
the \f$ x\f$ coordinate of `p` and the horizontal projection of `p` 
on `h`. 
*/ 
typedef unspecified_type Compare_x_at_y_2; 

/*!
Function object type that provides 
`Line_2 operator()(Point_2 p, Point_2 q)`, which constructs and 
returns the line defined by the points \f$ p\f$ and \f$ q\f$. 
*/ 
typedef unspecified_type Construct_line_2; 

/*!
Function object type that provides 
`bool operator()(Line_2 l)`, which returns `true` iff the 
line `l` is horizontal. 
*/ 
typedef unspecified_type Is_horizontal_2; 

/// @} 

/// \name Creation 
/// A copy constructor and default constructor are required.
/// @{

/*!

*/ 
YMonotonePartitionTraits(); 

/*!

*/ 
YMonotonePartitionTraits(const YMonotonePartitionTraits tr); 

/// @} 

/// \name Operations 
/// In addition to the functions required for the concept
/// `PartitionTraits_2`, the following functions that create instances
/// of the above function object types must exist.
/// @{

/*!

*/ 
Construct_line_2 construct_line_2_object(); 

/*!

*/ 
Compare_x_at_y_2 compare_x_at_y_2_object(); 

/*!

*/ 
Is_horizontal_2 is_horizontal_2_object(); 

/// @}

}; /* end YMonotonePartitionTraits_2 */
