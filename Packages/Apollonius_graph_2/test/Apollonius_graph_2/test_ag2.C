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

typedef CGAL::Simple_cartesian<number_t>   Kernel;

#include <CGAL/Number_type_traits.h>

typedef CGAL::Ring_tag Method_tag;

#include "./include/test.h"



int main()
{
  std::ifstream ifs_traits("./data/traits.dat");
  std::ifstream ifs_algo("./data/algo.dat");

  assert( ifs_traits );
  assert( ifs_algo );

  //  bool is_ok =
  //    CGAL::test_traits<Kernel,CGAL::Ring_tag,std::ifstream>(ifs_traits);

  bool is_ok = CGAL::test_traits<Kernel,Method_tag>();

  assert( is_ok );

  bool algo_ok =
    CGAL::test_algo<Kernel,Method_tag,std::ifstream>(ifs_algo);

  assert( algo_ok );

  return 0;
}
