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

typedef CGAL::Simple_cartesian<number_t>   Kernel;

#include <CGAL/Number_type_traits.h>

typedef CGAL::Ring_tag Method_tag;

#include "./include/test.h"



int main()
{
  std::ifstream ifs_traits("./data/traits.dat");
  std::ifstream ifs_algo("./data/algo.dat");
  std::ifstream ifs_hierarchy("./data/hierarchy.dat");

  assert( ifs_traits );
  assert( ifs_algo );
  assert( ifs_hierarchy );

  //  bool is_ok =
  //    CGAL::test_traits<Kernel,CGAL::Ring_tag,std::ifstream>(ifs_traits);

  bool traits_ok = CGAL::test_traits<Kernel,Method_tag>();

  assert( traits_ok );

  bool algo_ok =
    CGAL::test_algo<Kernel,Method_tag,std::ifstream>(ifs_algo);

  assert( algo_ok );

  bool hierarchy_ok = 
    CGAL::test_hierarchy_algo<Kernel,Method_tag,std::ifstream>(ifs_hierarchy);

  assert( hierarchy_ok );

  return 0;
}
