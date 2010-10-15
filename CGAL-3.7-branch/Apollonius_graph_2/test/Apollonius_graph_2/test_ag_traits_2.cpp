#include <CGAL/basic.h>

#include <iostream>
#include <fstream>
#include <cassert>

#define DONT_USE_FILTERED_EXACT

// choose number type
#include <CGAL/MP_Float.h>
#ifndef DONT_USE_FILTERED_EXACT
#  include <CGAL/Filtered_exact.h>
#endif

typedef double         inexact_type;
typedef CGAL::MP_Float exact_type;

#ifndef DONT_USE_FILTERED_EXACT
typedef CGAL::Filtered_exact<inexact_type,exact_type>  number_t;
#endif

#include <CGAL/Simple_cartesian.h>

#ifndef DONT_USE_FILTERED_EXACT
struct Kernel : public CGAL::Simple_cartesian<number_t> {};
#endif

typedef CGAL::Integral_domain_without_division_tag Method_tag;

#include "./include/test.h"


struct CK : public CGAL::Simple_cartesian<double> {};
struct EK : public CGAL::Simple_cartesian<CGAL::MP_Float> {};


int main()
{
#ifndef DONT_USE_FILTERED_EXACT
  {
    std::ifstream ifs_traits("./data/traits.dat");

    assert( ifs_traits );

    //  bool is_ok =
    //    CGAL::test_traits<Kernel,CGAL::Integral_domain_without_division_tag,std::ifstream>(ifs_traits);

    std::cout << "testing the traits class..." << std::flush;

    CGAL::Traits_tester<Kernel,Method_tag> test_traits;
    bool traits_ok = test_traits();

    assert( traits_ok );
    std::cout << " done!" << std::endl;

    ifs_traits.close();

    std::cout << std::endl;
  }
#endif
  //------------------------------------------------------------------------

  {
    std::ifstream ifs_traits("./data/traits.dat");

    assert( ifs_traits );

    std::cout << "testing the filtered traits class..." << std::flush;

    CGAL::Filtered_traits_tester<CK,Method_tag,EK,Method_tag> test_traits;
    bool traits_ok = test_traits();

    assert( traits_ok );
    std::cout << " done!" << std::endl;

    ifs_traits.close();

    std::cout << std::endl;
  }

  return 0;
}
