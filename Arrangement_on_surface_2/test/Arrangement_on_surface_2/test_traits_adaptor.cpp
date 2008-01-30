#include <CGAL/basic.h>
#include "test_configuration.h"
#include <iostream>

#include <CGAL/assertions.h>
#include <CGAL/Arrangement_2.h>

#include <vector>

#include "test_traits_adaptor.h"
#include "Traits_adaptor_test.h"


// Arrangement types:
typedef CGAL::Arr_default_dcel<Traits>                  Dcel;
typedef CGAL::Arrangement_2<Traits, Dcel>               Arr;

// Traits types:
typedef Traits::Point_2                                 Point_2;
typedef Traits::X_monotone_curve_2                      X_monotone_curve_2;
typedef Traits::Curve_2                                 Curve_2;

int main (int argc, char * argv[])
{
  CGAL::set_error_behaviour(CGAL::CONTINUE);
  CGAL::set_warning_behaviour(CGAL::CONTINUE);
  prev_error_handler = CGAL::set_error_handler(failure_handler);
  prev_warning_handler = CGAL::set_warning_handler(failure_handler);
  Traits_adaptor_test<Traits> test(argc, argv);
  bool rc;
  rc = test.start();
  CGAL::set_error_handler(prev_error_handler);
  CGAL::set_warning_handler(prev_warning_handler);
  return (rc) ? 0 : -1;
}

#if TEST_TRAITS == SPHERICAL_ARC_TRAITS

/*! Read a point */

template <>
template <class stream>
bool
Traits_adaptor_test<Traits >::
read_point(stream & is, Point_2 & p)
{
  Basic_number_type x, y, z;
  is >> x >> y >> z;
  p = Point_2(x, y, z);
  return true;
}

/*! Read a xcurve */
template <>
template <class stream>
bool
Traits_adaptor_test<Traits>::read_xcurve(stream & is, X_monotone_curve_2 & xcv)
{
  Point_2 p1,p2;
  read_point(is, p1);
  read_point(is, p2);
  CGAL_assertion(p1 != p2);
  xcv = X_monotone_curve_2(p1, p2);
  return true;
}

/*! Read a curve */
template <>
template <class stream>
bool
Traits_adaptor_test<Traits>::read_curve(stream & is, Curve_2 & cv)
{
  Point_2 p1, p2;
  read_point(is, p1);
  read_point(is, p2);
  CGAL_assertion(p1 != p2);
  cv = Curve_2(p1, p2);
  return true;
}

#endif
