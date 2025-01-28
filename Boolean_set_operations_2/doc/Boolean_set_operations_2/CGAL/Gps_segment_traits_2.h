
namespace CGAL {

/*!
\ingroup PkgBooleanSetOperations2Ref

\cgalModels{GeneralPolygonSetTraits_2}

\sa `CGAL::Arr_segment_traits_2<Kernel>`

*/
template< typename Kernel, typename Container, typename ArrSegmentTraits >
class Gps_segment_traits_2 {
public:

/// \name Definition
/// The traits class `Gps_segment_traits_2` models the concept
/// `GeneralPolygonSetTraits_2`. It enables Boolean set-operations on
/// (linear) polygons. It defines the exposed type `Polygon_2` to be
/// `CGAL::Polygon_2<Kernel,Container>`. By default, the template
/// parameter `Container` is instantiated by
/// `std::vector<Kernel::Point_2>` and the template parameter
/// `ArrSegmentTraits` is instantiated by
/// `CGAL::Arr_segment_traits_2<Kernel>`.
/// @{

/*!

*/
typedef CGAL::Polygon_2<Kernel,Container> Polygon_2;

/// @}

}; /* end Gps_segment_traits_2 */
} /* end namespace CGAL */
