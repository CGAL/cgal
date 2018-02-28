#ifndef TRAITS_LINESTRING_H
#define TRAITS_LINESTRING_H
#include <boost/geometry/geometry.hpp>
#include <CGAL/Geometry_container.h>



namespace boost{
namespace geometry{
namespace traits{
//!\todo should we use our own tag in namespace CGAL rather than use the ones from boost ?
template< typename R> struct tag<CGAL::Geometry_container<R, linestring_tag> >
{ typedef linestring_tag type; };

}}} //end namespaces

#endif // TRAITS_LINESTRING_H
