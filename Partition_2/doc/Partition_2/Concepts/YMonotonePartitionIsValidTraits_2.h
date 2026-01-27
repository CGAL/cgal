/*!
\ingroup PkgPartition2Concepts
\cgalConcept

Requirements of a traits class that is used
by `y_monotone_partition_is_valid_2` for testing the validity of a
\f$ y\f$-monotone partition of a polygon.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Partition_traits_2<R>}
\cgalHasModelsEnd

\sa `CGAL::partition_is_valid_2()`
\sa `CGAL::y_monotone_partition_2()`

*/

class YMonotonePartitionIsValidTraits_2 {
public:

/// \name Types
/// All types required by the concept `PartitionIsValidTraits_2` are
/// required except the function object type `Is_valid`. The following
/// type is required instead:
/// @{

/*!
Model of the concept `PolygonIsValid` that tests if
a sequence of points is \f$ y\f$-monotone or not.
*/
typedef unspecified_type Is_y_monotone_2;

/// @}

/// \name Creation
/// Only a copy constructor is required.
/// @{

/*!

*/
YMonotonePartitionIsValidTraits_2(YMonotonePartitionIsValidTraits_2 tr );

/// @}

/// \name Operations
/// The following function that creates an instance of the above
/// predicate object type must exist instead of the function
/// `is_valid_object` required by `PartitionIsValidTraits_2`.
/// @{

/*!

*/
Is_y_monotone_2 is_y_monotone_2_object(YMonotonePartitionIsValidTraits_2 tr);

/// @}

}; /* end YMonotonePartitionIsValidTraits_2 */
