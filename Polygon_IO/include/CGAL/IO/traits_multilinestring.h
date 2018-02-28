#ifndef TRAITS_MULTILINESTRING_H
#define TRAITS_MULTILINESTRING_H

#include <CGAL/Geometry_container.h>
#include <boost/geometry/geometry.hpp>

namespace boost{
namespace geometry{
namespace traits{
// WKT traits for MultiLinestring
template< typename R >
struct tag<CGAL::Geometry_container<R, multi_linestring_tag> >
{ typedef multi_linestring_tag type; };

}//end traits
}//end geometry
}//end boost
#endif // TRAITS_MULTILINESTRING_H
