#include "include/utils.h"
#include "include/wrappers.h"

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

using SCKER = CGAL::Simple_cartesian<double>;
using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;

template<typename Kernel>
void test_overloads()
{
  using FT      = typename Kernel::FT;
  using Point_2 = typename Kernel::Point_2;
  using Point_3 = typename Kernel::Point_3;

  const FT a = FT(1);
  const Point_2 p(0, 0);
  const Point_3 q(0, 0, 0);

  assert(CGAL::Weights::uniform_weight(p, p, p, p) == a);
  assert(CGAL::Weights::uniform_weight(q, q, q, q) == a);

  struct Traits : public Kernel { };
  assert(CGAL::Weights::uniform_weight(p, p, p, p, Traits()) == a);
  assert(CGAL::Weights::uniform_weight(q, q, q, q, Traits()) == a);
}

template<typename Kernel>
void test_kernel()
{
  test_overloads<Kernel>();
  const wrappers::Uniform_wrapper<Kernel> uni;
  tests::test_analytic_weight<Kernel>(uni, uni);
}

int main(int, char**)
{
  test_kernel<SCKER>();
  test_kernel<EPICK>();
  test_kernel<EPECK>();
  std::cout << "* test_uniform_weights: SUCCESS" << std::endl;
  return EXIT_SUCCESS;
}
