
namespace CGAL {

/*!
\ingroup PkgBooleanSetOperations2Ref

The class `Polygon_set_2` represents sets of linear polygons with holes.
The first two template parameters (`Kernel` and `Container`)
are used to instantiate the type `Polygon_2<Kernel,Container>`.
This type is used to represent the outer boundary of every set member
and the boundaries of all holes of every set member.

The third template parameter `Dcel` must be instantiated with a
model of the concept `GeneralPolygonSetDcel`. It is instantiated
by default with the type `Gps_default_dcel<Traits>`. You can override
this default, with a different \dcel class, typically an extension
of the `Gps_default_dcel` class template. Overriding the default is
necessary only if you intend to obtain the underlying internal arrangement
and process it further.

\sa `General_polygon_set_2`
\sa `Gps_segment_traits_2`

*/
template< typename Kernel, typename Container, typename Dcel >
class Polygon_set_2 : public General_polygon_set_2<Gps_segment_traits_2<Kernel,Container> > {
public:

}; /* end Polygon_set_2 */
} /* end namespace CGAL */
