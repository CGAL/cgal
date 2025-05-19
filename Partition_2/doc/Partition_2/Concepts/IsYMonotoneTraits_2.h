/*!
\ingroup PkgPartition2Concepts
\cgalConcept

Requirements of a traits class to be
used with the function `is_y_monotone_2()` that tests whether a sequence of
2D points defines a \f$ y\f$-monotone polygon or not.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Partition_traits_2<R>}
\cgalHasModels{CGAL::Kernel_traits_2}
\cgalHasModelsEnd

\sa `CGAL::Is_y_monotone_2<Traits>`
\sa `CGAL::y_monotone_partition_2()`
\sa `CGAL::y_monotone_partition_is_valid_2()`

*/

class IsYMonotoneTraits_2 {
public:

/// \name Types
/// The following two types are required:
/// @{

/*!
The point type of the polygon vertices.
*/
typedef unspecified_type Point_2;

/*!

Predicate object type that compares `Point_2`s lexicographically.
Must provide `bool operator()(Point_2 p, Point_2 q)` where `true`
is returned iff \f$ p <_{xy} q\f$.
We have \f$ p<_{xy}q\f$, iff \f$ p_x < q_x\f$ or \f$ p_x = q_x\f$ and \f$ p_y < q_y\f$,
where \f$ p_x\f$ and \f$ p_y\f$ denote \f$ x\f$ and \f$ y\f$ coordinate of point \f$ p\f$ resp.

*/
typedef unspecified_type Less_yx_2;

/// @}

/// \name Creation
/// Only a copy constructor is required.
/// @{

/*!

*/
IsYMonotoneTraits_2(IsYMonotoneTraits_2 tr);

/// @}

/// \name Operations
/// The following function that creates an instance of the above predicate object type must exist:
/// @{

/*!

*/
Less_yx_2 less_yx_2_object();

/// @}

}; /* end IsYMonotoneTraits_2 */
