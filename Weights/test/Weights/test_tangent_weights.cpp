#include "include/utils.h"
#include "include/wrappers.h"

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

using SCKER = CGAL::Simple_cartesian<double>;
using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;

template<typename Kernel>
void test_kernel()
{
  const wrappers::Tangent_wrapper<Kernel> tan;
  const wrappers::Mean_value_wrapper<Kernel> mvw;
  const wrappers::Three_point_family_wrapper<Kernel> tpf(1);
  tests::test_analytic_weight<Kernel>(tan, mvw);
  tests::test_analytic_weight<Kernel>(tan, tpf);
}

int main(int, char**)
{
  test_kernel<SCKER>();
  test_kernel<EPICK>();
  test_kernel<EPECK>();
  std::cout << "* test_tangent_weights: SUCCESS" << std::endl;
  return EXIT_SUCCESS;
}
