#ifndef CGAL_TRAITS_POLYGON_H
#define CGAL_TRAITS_POLYGON_H
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Point_2.h>
#include <boost/geometry/geometry.hpp>

namespace boost{
namespace geometry{
namespace traits{
// WKT traits for Polygon
template< typename K > struct tag<CGAL::Polygon_2<K> >
{ typedef ring_tag type; };

template< typename K >
struct tag<CGAL::Polygon_with_holes_2<K> >
{ typedef polygon_tag type; };

template< typename K >
struct ring_const_type<CGAL::Polygon_with_holes_2<K> >
{ typedef const CGAL::Polygon_2<K>& type; };

template< typename K >
struct ring_mutable_type<CGAL::Polygon_with_holes_2<K> >
{ typedef CGAL::Polygon_2<K>& type; };

template< typename K >
struct interior_const_type<CGAL::Polygon_with_holes_2<K> >
{ typedef const typename CGAL::Polygon_with_holes_2<K>::Holes_container& type; };

template< typename K >
struct interior_mutable_type<CGAL::Polygon_with_holes_2<K> >
{ typedef typename CGAL::Polygon_with_holes_2<K>::Holes_container& type; };

template< typename K >
struct exterior_ring<CGAL::Polygon_with_holes_2<K> >
{
  static CGAL::Polygon_2<K>& get(CGAL::Polygon_with_holes_2<K>& p)
  {
    return (p.outer_boundary());
  }
  static CGAL::Polygon_2<K> const& get(CGAL::Polygon_with_holes_2<K> const& p)
  {
    return (p.outer_boundary());
  }
};

template< typename K >
struct interior_rings<CGAL::Polygon_with_holes_2<K> >
{
  static typename CGAL::Polygon_with_holes_2<K>::Holes_container& get(CGAL::Polygon_with_holes_2<K>& p)
  {
    return p.holes();
  }
  static const typename CGAL::Polygon_with_holes_2<K>::Holes_container& get(CGAL::Polygon_with_holes_2<K> const& p)
  {
    return p.holes();
  }
};
}//end traits
}//end geometry

//extra specialization
template< typename K >
struct range_value<CGAL::Polygon_2<K> >
{
  typedef typename CGAL::Polygon_2<K>::Point_2  type;
};

}//end boost

#endif 