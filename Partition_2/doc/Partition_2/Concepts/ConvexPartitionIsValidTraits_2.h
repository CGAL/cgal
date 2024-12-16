/*!
\ingroup PkgPartition2Concepts
\cgalConcept

Requirements of a traits class used
by `convex_partition_is_valid_2` for testing the validity of a
convex partition of a polygon.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Partition_traits_2<R>}
\cgalHasModelsEnd

\sa `CGAL::approx_convex_partition_2()`
\sa `CGAL::greene_approx_convex_partition_2()`
\sa `CGAL::Is_convex_2<Traits>`
\sa `CGAL::optimal_convex_partition_2()`

*/

class ConvexPartitionIsValidTraits_2 {
public:

/// \name Types
/// All types required by the concept `PartitionIsValidTraits_2` are
/// required except the function object type `Is_valid`. The following
/// type is required instead: @{
/*!
Model of the concept `PolygonIsValid` that tests if
a sequence of points is convex or not.
*/
typedef unspecified_type Is_convex_2;

/// @}

/// \name Creation
/// Only a copy constructor is required.
/// @{

/*!

*/
ConvexPartitionIsValidTraits_2(ConvexPartitionIsValidTraits_2 tr);

/// @}

/// \name Operations
/// The following function that creates an instance of the above
/// predicate object type must exist instead of the function
/// `is_valid_object` required by `PartitionIsValidTraits_2`.
/// @{

/*!

*/
Is_convex_2 is_convex_2_object(ConvexPartitionIsValidTraits_2 tr);

/// @}

}; /* end ConvexPartitionIsValidTraits_2 */
