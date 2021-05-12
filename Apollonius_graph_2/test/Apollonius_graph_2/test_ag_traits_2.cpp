#include <iostream>
#include <fstream>
#include <cassert>

// choose number type
#include <CGAL/MP_Float.h>

typedef double         inexact_type;
typedef CGAL::MP_Float exact_type;

#include <CGAL/Simple_cartesian.h>

typedef CGAL::Integral_domain_without_division_tag Method_tag;

#include "./include/test.h"

typedef CGAL::Simple_cartesian<double> CK;
typedef CGAL::Simple_cartesian<CGAL::MP_Float> EK;

int main()
{
  std::ifstream ifs_traits("./data/traits.dat");
  assert( ifs_traits );

  std::cout << "testing the filtered traits class..." << std::flush;

  CGAL::Filtered_traits_tester<CK,Method_tag,EK,Method_tag> test_traits;
  bool traits_ok = test_traits();
  assert( traits_ok );

  std::cout << " done!" << std::endl;

  return 0;
}
