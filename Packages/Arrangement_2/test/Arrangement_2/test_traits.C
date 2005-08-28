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
bool
Traits_test<CGAL::Arr_segment_traits_2<Kernel> >::
read_curve(std::ifstream & is,
           CGAL::Arr_segment_traits_2<Kernel>::Curve_2 & cv)
{
  return false;
}

template <>
bool
Traits_test<CGAL::Arr_segment_traits_2<Kernel> >::
read_point(std::ifstream & is,
           CGAL::Arr_segment_traits_2<Kernel>::Point_2 & p)
{
  Number_type x, y;
  is >> x >> y;
  p = Point_2(x, y);
  return true;
}
