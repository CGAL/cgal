#include <CGAL/basic.h>

// Workaround for buggy compilers.
#ifdef CGAL_CFG_MATCHING_BUG_2
#  define CGAL_IA_CT double
#  define CGAL_IA_PROTECTED true
#  define CGAL_IA_CACHE No_Filter_Cache
#  define CGAL_IA_ET CGAL::MP_Float
#endif

#include <iostream>
#include <fstream>
#include <cassert>

// choose number type
#include <CGAL/MP_Float.h>
#include <CGAL/Filtered_exact.h>

typedef double         inexact_type;
typedef CGAL::MP_Float exact_type;

typedef CGAL::Filtered_exact<inexact_type,exact_type>  number_t;

#include <CGAL/Simple_cartesian.h>

struct Kernel : public CGAL::Simple_cartesian<number_t> {};

#include <CGAL/Number_type_traits.h>

typedef CGAL::Ring_tag Method_tag;

#include "./include/test.h"


struct CK : public CGAL::Simple_cartesian<double> {};
struct EK : public CGAL::Simple_cartesian<CGAL::MP_Float> {};


int main()
{
  {
    std::ifstream ifs_traits("./data/traits.dat");
    std::ifstream ifs_algo("./data/algo.dat");
    std::ifstream ifs_hierarchy("./data/hierarchy.dat");

    assert( ifs_traits );
    assert( ifs_algo );
    assert( ifs_hierarchy );

    //  bool is_ok =
    //    CGAL::test_traits<Kernel,CGAL::Ring_tag,std::ifstream>(ifs_traits);

    std::cout << "testing the traits class..." << std::flush;

    CGAL::Traits_tester<Kernel,Method_tag> test_traits;
    bool traits_ok = test_traits();

    assert( traits_ok );
    std::cout << " done!" << std::endl;

    std::cout << "testing the Apollonius graph class..." << std::flush;
    bool algo_ok =
      CGAL::test_algo<Kernel,Method_tag,std::ifstream>(ifs_algo);

    assert( algo_ok );
    std::cout << " done!" << std::endl;

    std::cout << "testing the Apollonius graph hierarchy class..."
	      << std::flush;
    bool hierarchy_ok = 
      CGAL::test_hierarchy_algo<Kernel,Method_tag,
      std::ifstream>(ifs_hierarchy);

    assert( hierarchy_ok );
    std::cout << " done!" << std::endl;

    ifs_traits.close();
    ifs_algo.close();
    ifs_hierarchy.close();

    std::cout << std::endl;
  }

  //------------------------------------------------------------------------

  {
    std::ifstream ifs_traits("./data/traits.dat");
    std::ifstream ifs_algo("./data/algo.dat");
    std::ifstream ifs_hierarchy("./data/hierarchy.dat");

    assert( ifs_traits );
    assert( ifs_algo );
    assert( ifs_hierarchy );

    std::cout << "testing the filtered traits class..." << std::flush;

    CGAL::Filtered_traits_tester<CK,Method_tag,EK,Method_tag> test_traits;
    bool traits_ok = test_traits();

    assert( traits_ok );
    std::cout << " done!" << std::endl;

    std::cout << "testing the Apollonius graph class"
	      << " with filtered traits..." << std::flush;
    bool algo_ok =
      CGAL::test_filtered_traits_algo<CK,Method_tag,EK,Method_tag,
      std::ifstream>(ifs_algo);

    assert( algo_ok );
    std::cout << " done!" << std::endl;

    std::cout << "testing the Apollonius graph hierarchy class"
	      << " with filtered traits..." << std::flush;
    bool hierarchy_ok = 
      CGAL::test_filtered_traits_hierarchy_algo<CK,Method_tag,EK,
      Method_tag,std::ifstream>(ifs_hierarchy);

    assert( hierarchy_ok );
    std::cout << " done!" << std::endl;

    ifs_traits.close();
    ifs_algo.close();
    ifs_hierarchy.close();

    std::cout << std::endl;
  }

  return 0;
}
