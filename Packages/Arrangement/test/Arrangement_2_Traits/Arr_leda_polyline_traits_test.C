#include <CGAL/basic.h>
#include <iostream>

// Making sure test doesn't fail if LEDA is not installed
#if ! defined(CGAL_USE_LEDA)

int main()
{
  std::cout << "A try to run test with LEDA traits but LEDA is not installed.";
  std::cout << std::endl;
  std::cout << "Test is not performed.";
  std::cout << std::endl;

  return 0;
}
#else

#include <CGAL/leda_rational.h>
#include <CGAL/Arr_leda_polyline_traits.h>
#include <CGAL/Pm_segment_traits_leda_kernel_2.h>
#include "include/Polyline_traits_test.h"

typedef leda_rational                                   NT;
typedef CGAL::Pm_segment_traits_leda_kernel_2           Kernel;
typedef CGAL::Arr_leda_polyline_traits<Kernel>          Traits;

int main(int argc, char * argv[])
{
  Polyline_traits_test<Traits, leda_rational> test_obj(argc, argv);
  return (test_obj.start()) ? 0 /* SUCCESS */ : 1 /* FAILURE */;
}

#endif
