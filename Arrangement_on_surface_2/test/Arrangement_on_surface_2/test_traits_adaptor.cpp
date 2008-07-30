#include <CGAL/basic.h>
#include "test_configuration.h"
#include <iostream>

#include <CGAL/assertions.h>
#include <CGAL/Arrangement_2.h>

#include "test_traits_adaptor.h"
#include "Traits_adaptor_test.h"

int main (int argc, char * argv[])
{
  CGAL::set_error_behaviour(CGAL::THROW_EXCEPTION);
  CGAL::set_warning_behaviour(CGAL::THROW_EXCEPTION);
  Traits_adaptor_test<Traits> test(argc, argv);
  bool rc = test.start();
  return (rc) ? 0 : -1;
}
