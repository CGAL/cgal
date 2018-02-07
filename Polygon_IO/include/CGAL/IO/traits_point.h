#ifndef CGAL_TRAITS_POINT_H
#define CGAL_TRAITS_POINT_H

#include <CGAL/Point_2.h>
#include <boost/geometry/geometry.hpp>
namespace boost{
namespace geometry{
namespace traits{

//Traits for Points
template< typename K > struct tag<CGAL::Point_2<K> >
{ typedef point_tag type; };

template< typename K > struct coordinate_type<CGAL::Point_2<K> >
{ typedef typename K::FT type; };

template< typename K > struct coordinate_system<CGAL::Point_2<K> >
{ typedef cs::cartesian type; };

template< typename K > struct dimension<CGAL::Point_2<K> > : boost::mpl::int_<2> {};

template< typename K >
struct access<CGAL::Point_2<K> , 0>
{
  static double get(CGAL::Point_2<K>  const& p)
  {
    return CGAL::to_double(p.x());
  }
  
  static void set(CGAL::Point_2<K> & p, typename K::FT c)
  {
    p = CGAL::Point_2<K> (c, p.y());
  }
  
};

template< typename K >
struct access<CGAL::Point_2<K> , 1>
{
  static double get(CGAL::Point_2<K>  const& p)
  {
    return CGAL::to_double(p.y());
  }
  
  static void set(CGAL::Point_2<K> & p, typename K::FT c)
  {
    p = CGAL::Point_2<K> (p.x(), c);
  }
  
};

}}}//end namespaces
#endif