#ifndef CGAL_TRAITS_MULTIPOLYGON_H
#define CGAL_TRAITS_MULTIPOLYGON_H
#include <CGAL/Geometry_container.h>
#include <boost/geometry/geometry.hpp>


namespace boost{
namespace geometry{
namespace traits{
// WKT traits for MultiPolygon
template< typename R >
struct tag<CGAL::Geometry_container<R, multi_polygon_tag> >
{ typedef multi_polygon_tag type; };

}//end traits
}//end geometry
}//end boost

#endif 

