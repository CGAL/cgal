#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>

// Making sure test doesn't fail if CORE is not installed
#ifndef CGAL_USE_CORE

int main(int argc, char* argv[])
{
  std::cout << "A try to run test with LEDA traits but LEDA is not installed.";
  std::cout << std::endl;
  std::cout << "Test is not performed.";
  std::cout << std::endl;

  return 0;
}
#else

#include <CORE/BigInt.h>
#include <CGAL/CORE_Expr.h>
#include <CGAL/Arr_2_bases.h>
#include <CGAL/Arr_2_default_dcel.h>
#include <CGAL/Arrangement_2.h>

#include <CGAL/Arr_conic_traits_2.h>
#include "include/Conic_traits_test.h"

typedef CORE::BigInt                                  CfNT;
typedef CGAL::Cartesian<CfNT>                         IKernel;
typedef CORE::Expr                                    CoNT;
typedef CGAL::Cartesian<CoNT>                         AKernel;
typedef CGAL::Arr_conic_traits_2<IKernel, AKernel>    Traits;

int main (int argc, char** argv)
{
  Conic_traits_test<Traits, CoNT>  test_obj (argc, argv);

  test_obj.start();
  return (0);
}

#endif
