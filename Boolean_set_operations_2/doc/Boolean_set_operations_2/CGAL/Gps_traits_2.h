
namespace CGAL {

/*!
\ingroup PkgBooleanSetOperations2Ref

The traits class `Gps_traits_2` models the concept
`GeneralPolygonSetTraits_2`. It inherits from the instantiated
type of the template parameter `ArrTraits`, which must model the
concept `ArrangementDirectionalXMonotoneTraits_2`, (which in turn refines
the concept `ArrangementXMonotoneTraits_2`). The template parameter
`GeneralPolygon_t` must be instantiated with a model of the concept
of `GpsTraitsGeneralPolygon_2`. By default, the latter is instantiated by
`CGAL::General_polygon_2<ArrTraits>`.

\cgalModels{GeneralPolygonSetTraits_2}

*/
template< typename ArrTraits, typename GeneralPolygon_t >
class Gps_traits_2 {
public:

}; /* end Gps_traits_2 */
} /* end namespace CGAL */
