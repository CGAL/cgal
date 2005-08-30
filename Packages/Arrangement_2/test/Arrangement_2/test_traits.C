#include <CGAL/basic.h>
#include <CGAL/Arrangement_2.h>

#include "test_configuration.h"
#include "test_traits.h"
#include "Traits_test.h"

// Arrangement types:
typedef CGAL::Arr_default_dcel<Traits>                  Dcel;
typedef CGAL::Arrangement_2<Traits, Dcel>               Arr;

// Traits types:
typedef Traits::Point_2                                 Point_2;
typedef Traits::X_monotone_curve_2                      X_monotone_curve_2;
typedef Traits::Curve_2                                 Curve_2;

int main (int argc, char * argv[])
{
  Traits_test<Traits> test(argc, argv);
  return (test.start()) ? 0 /* SUCCESS */ : -1; /* FAILURE */
}

template <>
template <class stream>
bool
Traits_test<CGAL::Arr_segment_traits_2<Kernel> >::
read_point(stream & is,
           CGAL::Arr_segment_traits_2<Kernel>::Point_2 & p)
{
  Basic_number_type x, y;
  is >> x >> y;
  p = CGAL::Arr_segment_traits_2<Kernel>::Point_2(x, y);
  return true;
}

template <>
template <class stream>
bool
Traits_test<CGAL::Arr_segment_traits_2<Kernel> >::
read_xcurve(stream & is,
            CGAL::Arr_segment_traits_2<Kernel>::X_monotone_curve_2 & cv)
{
  Basic_number_type x1, y1, x2, y2;
  is >> x1 >> y1 >> x2 >> y2;
  Point_2 p1(x1, y1);
  Point_2 p2(x2, y2);
  CGAL_assertion(p1 != p2);
  cv = CGAL::Arr_segment_traits_2<Kernel>::X_monotone_curve_2(p1, p2);
  return true;
}

template <>
template <class stream>
bool
Traits_test<CGAL::Arr_segment_traits_2<Kernel> >::
read_curve(stream & is,
           CGAL::Arr_segment_traits_2<Kernel>::Curve_2 & cv)
{
  Basic_number_type x1, y1, x2, y2;
  is >> x1 >> y1 >> x2 >> y2;
  Point_2 p1(x1, y1);
  Point_2 p2(x2, y2);
  CGAL_assertion(p1 != p2);
  cv = CGAL::Arr_segment_traits_2<Kernel>::Curve_2(p1, p2);
  return true;
}

