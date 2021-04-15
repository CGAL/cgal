namespace CGAL {

/*!
\ingroup PkgPartition2Ref

Traits class that can be used with all the
2-dimensional polygon partitioning algorithms.


\tparam R  a representation class
\tparam PointPropertyMap a property map that maps to points of type `R::Point_2`

\cgalModels `ConvexPartitionIsValidTraits_2`
\cgalModels `IsYMonotoneTraits_2`
\cgalModels `OptimalConvexPartitionTraits_2`
\cgalModels `PartitionTraits_2`
\cgalModels `YMonotonePartitionIsValidTraits_2`

\sa `CGAL::approx_convex_partition_2()`
\sa `CGAL::convex_partition_is_valid_2()`
\sa `CGAL::greene_approx_convex_partition_2()`
\sa `CGAL::optimal_convex_partition_2()`
\sa `CGAL::partition_is_valid_2()`
\sa `CGAL::Partition_is_valid_traits_2<Traits, PolygonIsValid>`
\sa `CGAL::y_monotone_partition_2()`
\sa `CGAL::y_monotone_partition_is_valid_2()`

*/
  template< typename R, typename PointPropertyMap = Identity_property_map<R::Point_2> >
class Partition_traits_2 {
public:

/// \name Types
/// @{

/*!

*/
    typedef Partition_traits_2<R,PMap> Self;

/*!

*/
typedef R::FT FT;

/*!

*/
typedef boost::property_traits<PointPropertyMap>::key_type Point_2;


/*!

*/
typedef std::list<Point_2> Container;

/*!

*/
typedef CGAL::Polygon_2<Self, Container> Polygon_2;



/*!
A functor with an operator which first obtains points of type `R::Point_2`
with the function `get()` applied on the point property map, and
then applies the functor of `R::Less_yx_2` to these points.
*/
typedef unspecified_type Less_yx_2;

/*!

*/
typedef unspecified_type Less_xy_2;

/*!

*/
typedef unspecified_type Left_turn_2;

/*!

*/
typedef unspecified_type Orientation_2;

/*!

*/
typedef unspecified_type Compare_y_2;

/*!

*/
typedef unspecified_type Compare_x_2;


/*!

*/
typedef unspecified_type Collinear_are_ordered_along_line_2;

/*!

*/
    typedef unspecified_type Are_strictly_ordered_along_line_2;

/*!
*/
typedef CGAL::Is_convex_2<Self> Is_convex_2;

/*!

*/
typedef CGAL::Is_y_monotone_2<Self> Is_y_monotone_2;

/// @}

/// \name Creation
/// A default constructor and copy constructor are defined.
/// @{

/*!

*/
Partition_traits_2();

/*!

*/
Partition_traits_2(Partition_traits_2& tr);

/*!

*/
Partition_traits_2(const R& rep);

/*!

*/
Partition_traits_2(const R& rep, PointPropertyMap pmap);

/// @}

/// \name Operations
/// For each predicate object type `Pred_object_type` listed above
/// (i.e., `Less_yx_2`, `Less_xy_2`, `Left_turn_2`,
/// `Orientation_2`, `Compare_y_2`, `Compare_x_2`,
/// `Collinear_are_ordered_along_line_2`,
/// `Are_strictly_ordered_along_line_2`, `Is_convex_2`,
/// `Is_y_monotone_2`) there is a corresponding function of the
/// following form defined:
/// @{

/*!

Returns an instance of `Pred_object_type`.

*/
Pred_object_type pred_object_type_object();

/// @}

}; /* end Partition_traits_2 */
} /* end namespace CGAL */
