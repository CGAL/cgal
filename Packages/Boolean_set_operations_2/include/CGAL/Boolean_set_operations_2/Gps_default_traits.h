
#ifndef GPS_DEFAULT_TRAITS_H
#define GPS_DEFAULT_TRAITS_H

#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/General_polygon_with_holes_2.h>
#include <CGAL/General_polygon_2.h>
#include <CGAL/Gps_segment_traits_2.h>
#include <CGAL/Gps_traits_2.h>

CGAL_BEGIN_NAMESPACE

template <class Polygon>
struct Gps_default_traits
{};


template <class Kernel, class Container>
struct Gps_default_traits<CGAL::Polygon_2<Kernel, Container> >
{
  typedef Gps_segment_traits_2<Kernel,
                               Container,
                               Arr_segment_traits_2<Kernel> >    Traits;
};

template <class Kernel, class Container>
struct Gps_default_traits<CGAL::Polygon_with_holes_2<Kernel, Container> >
{
  typedef Gps_segment_traits_2<Kernel,
                               Container,
                               Arr_segment_traits_2<Kernel> >    Traits;
};

template <class Polygon>
struct Gps_default_traits<CGAL::General_polygon_with_holes_2<Polygon> >
{
  typedef typename Gps_default_traits<Polygon>::Traits Traits;
};

template <class Arr_traits>
struct Gps_default_traits<CGAL::General_polygon_2<Arr_traits> >
{
  typedef Gps_traits_2<Arr_traits>    Traits;
};

CGAL_END_NAMESPACE

#endif
