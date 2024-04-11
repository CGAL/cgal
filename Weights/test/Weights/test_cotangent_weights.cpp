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
  using FT      = typename Kernel::FT;
  using Point_2 = typename Kernel::Point_2;

  const Point_2 p0(-2, 1);
  const Point_2 p1( 0, 1);
  const Point_2 p2( 0, 3);
  const Point_2 q( -2, 3);
  const FT w = CGAL::Weights::cotangent_weight(p0, p1, p2, q);
  assert(w == FT(0));

  const wrappers::Cotangent_wrapper<Kernel> cot;
  const wrappers::Discrete_harmonic_wrapper<Kernel> dhw;
  const wrappers::Three_point_family_wrapper<Kernel> tpf(2);
  tests::test_analytic_weight<Kernel>(cot, dhw);
  tests::test_analytic_weight<Kernel>(cot, tpf);
}

int main(int, char**)
{
  test_kernel<SCKER>();
  test_kernel<EPICK>();
  test_kernel<EPECK>();
  std::cout << "* test_cotangent_weights: SUCCESS" << std::endl;
  return EXIT_SUCCESS;
}
