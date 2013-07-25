namespace CGAL {

/*!
\ingroup PkgPolygonPartitioning2

Traits class that can be used with all the 
2-dimensional polygon partitioning algorithms. It is parameterized by 
a representation class `R`. 

\cgalModels `ConvexPartitionIsValidTraits_2`
\cgalModels `IsYMonotoneTraits_2`
\cgalModels `OptimalConvexPartitionTraits_2`
\cgalModels `PartitionTraits_2`
\cgalModels `YMonotonePartitionIsValidTraits_2`
\cgalModels `YMonotonePartitionTraits_2`

\sa `CGAL::approx_convex_partition_2()` 
\sa `CGAL::convex_partition_is_valid_2()` 
\sa `CGAL::greene_approx_convex_partition_2()` 
\sa `CGAL::optimal_convex_partition_2()` 
\sa `CGAL::partition_is_valid_2()` 
\sa `CGAL::Partition_is_valid_traits_2<Traits, PolygonIsValid>` 
\sa `CGAL::y_monotone_partition_2()` 
\sa `CGAL::y_monotone_partition_is_valid_2()` 

*/
template< typename R >
class Partition_traits_2 {
public:

/// \name Types 
/// @{

/*!

*/ 
typedef R::Line_2 Line_2; 

/*!

*/ 
typedef R::Segment_2 Segment_2; 

/*!

*/ 
typedef R::Ray_2 Ray_2; 

/*!

*/ 
typedef R::Less_yx_2 Less_yx_2; 

/*!

*/ 
typedef R::Less_xy_2 Less_xy_2; 

/*!

*/ 
typedef R::Left_turn_2 Left_turn_2; 

/*!

*/ 
typedef R::Orientation_2 Orientation_2; 

/*!

*/ 
typedef R::Compare_y_2 Compare_y_2; 

/*!

*/ 
typedef R::Compare_x_2 Compare_x_2; 

/*!

*/ 
typedef R::Construct_line_2 Construct_line_2; 

/*!

*/ 
typedef R::Construct_ray_2 Construct_ray_2; 

/*!

*/ 
typedef R::Construct_segment_2 Construct_segment_2; 

/*!

*/ 
typedef R::Collinear_are_ordered_along_line_2 Collinear_are_ordered_along_line_2; 

/*!

*/ 
typedef R::Are_strictly_ordered_along_line_2 Are_strictly_ordered_along_line_2; 

/*!

*/ 
typedef CGAL::Polygon_traits_2<R> Poly_Traits; 

/*!

*/ 
typedef Poly_Traits::Point_2 Point_2; 

/*!

*/ 
typedef std::list<Point_2> Container; 

/*!

*/ 
typedef CGAL::Polygon_2<Poly_Traits, Container> Polygon_2; 

/*!

*/ 
typedef R::Less_xy_2 Less_xy; 

/*!

*/ 
typedef Poly_Traits::Vector_2 Vector_2; 

/*!

*/ 
typedef R::FT FT; 

/*!

*/ 
typedef Partition_traits_2<R> Self; 

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

/// @} 

/// \name Operations 
/// For each predicate object type `Pred_object_type` listed above
/// (i.e., `Less_yx_2`, `Less_xy_2`, `Left_turn_2`,
/// `Orientation_2`, `Compare_y_2`, `Compare_x_2`, `Construct_line_2`,
/// `Construct_ray_2`, `Construct_segment_2`,
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
