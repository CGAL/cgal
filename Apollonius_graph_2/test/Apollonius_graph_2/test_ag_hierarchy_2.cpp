#include <CGAL/basic.h>

#include <iostream>
#include <fstream>
#include <cassert>

#define DNOT_USE_FILTERED_EXACT

// choose number type
#include <CGAL/MP_Float.h>

#ifndef DNOT_USE_FILTERED_EXACT
#  include <CGAL/Filtered_exact.h>
#endif

typedef double         inexact_type;
typedef CGAL::MP_Float exact_type;

#ifndef DNOT_USE_FILTERED_EXACT
typedef CGAL::Filtered_exact<inexact_type,exact_type>  number_t;
#endif

#include <CGAL/Simple_cartesian.h>

#ifndef DNOT_USE_FILTERED_EXACT
struct Kernel : public CGAL::Simple_cartesian<number_t> {};
#endif

typedef CGAL::Integral_domain_without_division_tag Method_tag;

#include "./include/test.h"


struct CK : public CGAL::Simple_cartesian<inexact_type> {};
struct EK : public CGAL::Simple_cartesian<exact_type> {};


int main()
{
#ifndef DNOT_USE_FILTERED_EXACT
  {
    std::ifstream ifs_hierarchy("./data/hierarchy.dat");

    assert( ifs_hierarchy );

    std::cout << "testing the Apollonius graph hierarchy class..."
	      << std::flush;
    bool hierarchy_ok = 
      CGAL::test_hierarchy_algo<Kernel,Method_tag,
      std::ifstream>(ifs_hierarchy);

    assert( hierarchy_ok );
    std::cout << " done!" << std::endl;

    ifs_hierarchy.close();

    std::cout << std::endl;
  }
#endif
  //------------------------------------------------------------------------

  {
    std::ifstream ifs_hierarchy("./data/hierarchy.dat");

    assert( ifs_hierarchy );

    std::cout << "testing the Apollonius graph hierarchy class"
	      << " with filtered traits..." << std::flush;
    bool hierarchy_ok = 
      CGAL::test_filtered_traits_hierarchy_algo<CK,Method_tag,EK,
      Method_tag,std::ifstream>(ifs_hierarchy);

    assert( hierarchy_ok );
    std::cout << " done!" << std::endl;

    ifs_hierarchy.close();

    std::cout << std::endl;
  }

  return 0;
}
