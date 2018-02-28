#ifndef TRAITS_MULTIPOINT_H
#define TRAITS_MULTIPOINT_H
#include <boost/geometry/geometry.hpp>
#include <CGAL/Geometry_container.h>

namespace boost{
namespace geometry{
namespace traits{
// WKT traits for MultiPoint
template< typename R >
struct tag<CGAL::Geometry_container<R, multi_point_tag > >
{ typedef multi_point_tag type; };

}//end traits
}//end geometry
}//end boost
#endif // TRAITS_MULTIPOINT_H
