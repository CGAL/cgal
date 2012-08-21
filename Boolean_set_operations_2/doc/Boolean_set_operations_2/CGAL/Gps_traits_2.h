
namespace CGAL {

/*!
\ingroup PkgBooleanSetOperations2

The traits class `Gps_traits_2` models the concept 
`GeneralPolygonSetTraits_2`. It inherits from the instantiated 
type of the template parameter `ArrTraits`, which must model the 
concept `ArrangementDirectionalXMonotoneTraits`, (which in turn refines 
the concept `ArrangementXMonotoneTraits`). The template parameter 
`GeneralPolygon_t` must be instantiated with a model of the concept 
of `GpsTraitsGeneralPolygon_2`. By default, the latter is instantiated by 
`CGAL::General_polygon_2<ArrTraits>`. 

\models ::GeneralPolygonSetTraits_2 

*/
template< typename ArrTraits, typename GeneralPolygon_t >
class Gps_traits_2 {
public:

/// @}

}; /* end Gps_traits_2 */
} /* end namespace CGAL */
