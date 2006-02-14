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
    std::ifstream ifs_algo("./data/algo.dat");

    assert( ifs_algo );

    std::cout << "testing the Apollonius graph class..." << std::flush;
    bool algo_ok =
      CGAL::test_algo<Kernel,Method_tag,std::ifstream>(ifs_algo);

    assert( algo_ok );
    std::cout << " done!" << std::endl;

    ifs_algo.close();

    std::cout << std::endl;
  }

  //------------------------------------------------------------------------

  {
    std::ifstream ifs_algo("./data/algo.dat");

    assert( ifs_algo );

    std::cout << "testing the Apollonius graph class"
	      << " with filtered traits..." << std::flush;
    bool algo_ok =
      CGAL::test_filtered_traits_algo<CK,Method_tag,EK,Method_tag,
      std::ifstream>(ifs_algo);

    assert( algo_ok );
    std::cout << " done!" << std::endl;

    ifs_algo.close();

    std::cout << std::endl;
  }

  return 0;
}
