#include <CGAL/basic.h>
#include "test_configuration.h"
#include <iostream>

#include <CGAL/assertions.h>
#include <CGAL/Arrangement_2.h>

#include <vector>

#include "test_traits.h"
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
