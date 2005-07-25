#include <CGAL/basic.h>

#include <iostream>
#include <fstream>
#include <cassert>

// choose number type
#ifdef CGAL_USE_GMP
#  include <CGAL/Gmpq.h>
#else
#  include <CGAL/MP_Float.h>
#endif

#include <CGAL/Filtered_exact.h>

typedef double         inexact_type;
#ifdef CGAL_USE_GMP
typedef CGAL::Gmpq     exact_type;
#else
typedef CGAL::MP_Float exact_type;
#endif

typedef CGAL::Filtered_exact<inexact_type,exact_type>  number_t;

#include <CGAL/Simple_cartesian.h>

struct Kernel : public CGAL::Simple_cartesian<number_t> {};

#include <CGAL/Number_type_traits.h>

typedef CGAL::Ring_tag Method_tag;

#include "./include/test.h"


struct CK : public CGAL::Simple_cartesian<double> {};
struct EK : public CGAL::Simple_cartesian<exact_type> {};


int main(int argc, char* argv[])
{
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
