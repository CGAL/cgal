
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

typedef CGAL::Simple_cartesian<inexact_type> CK;
typedef CGAL::Simple_cartesian<exact_type> EK;

int main()
{
  std::ifstream ifs_hierarchy("./data/hierarchy.dat");
  assert( ifs_hierarchy );

  std::cout << "testing the Apollonius graph hierarchy class" << " with filtered traits..." << std::flush;
  bool hierarchy_ok = CGAL::test_filtered_traits_hierarchy_algo<CK,Method_tag,EK,Method_tag,std::ifstream>(ifs_hierarchy);

  assert( hierarchy_ok );
  std::cout << " done!" << std::endl;

  return 0;
}
