#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>

// Making sure test doesn't fail if LEDA is not installed
#ifndef CGAL_USE_LEDA

int main(int argc, char* argv[])
{
  std::cout << "A try to run test with LEDA traits but LEDA is not installed.";
  std::cout << std::endl;
  std::cout << "Test is not performed.";
  std::cout << std::endl;

  return 0;
}
#else

#include <CGAL/leda_real.h>
#include <CGAL/Arr_2_bases.h>
#include <CGAL/Arr_2_default_dcel.h>
#include <CGAL/Arrangement_2.h>

#include <CGAL/Arr_conic_traits_2.h>
#include "include/Conic_traits_test.h"

typedef leda_real                         NT;
typedef CGAL::Arr_conic_traits_2<NT>      Traits;

int main (int argc, char** argv)
{
  Conic_traits_test<Traits, NT>  test_obj (argc, argv);

  test_obj.start();
  return (0);
}

#endif
