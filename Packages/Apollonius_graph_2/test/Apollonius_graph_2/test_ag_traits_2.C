#include <CGAL/basic.h>

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


int main(int argc, char* argv[])
{
  {
    std::ifstream ifs_traits("./data/traits.dat");

    assert( ifs_traits );

    //  bool is_ok =
    //    CGAL::test_traits<Kernel,CGAL::Ring_tag,std::ifstream>(ifs_traits);

    std::cout << "testing the traits class..." << std::flush;

    CGAL::Traits_tester<Kernel,Method_tag> test_traits;
    bool traits_ok = test_traits();

    assert( traits_ok );
    std::cout << " done!" << std::endl;

    ifs_traits.close();

    std::cout << std::endl;
  }

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
