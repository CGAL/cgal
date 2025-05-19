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

  const Point_2 p1(0, 0);
  const Point_2 q1(2, 0);
  const Point_3 p2(0, 0, 1);
  const Point_3 q2(2, 0, 1);

  const FT a2 = CGAL::Weights::shepard_weight(p1, q1, 3);
  const FT a3 = CGAL::Weights::shepard_weight(p2, q2, 3);
  assert(a2 == FT(1)/FT(8));
  assert(a3 == FT(1)/FT(8));

  assert(CGAL::Weights::shepard_weight(p1, p1, q1, q1, 3) == a2);
  assert(CGAL::Weights::shepard_weight(p2, p2, q2, q2, 3) == a3);

  struct Traits : public Kernel { };
  assert(CGAL::Weights::shepard_weight(p1, p1, q1, q1, 3, Traits()) == a2);
  assert(CGAL::Weights::shepard_weight(p2, p2, q2, q2, 3, Traits()) == a3);
}

template<typename Kernel>
void test_kernel()
{
  test_overloads<Kernel>();
  const wrappers::Shepard_wrapper<Kernel> spwa(1);
  const wrappers::Shepard_wrapper<Kernel> spwb(2);
  const wrappers::Inverse_distance_wrapper<Kernel> idw;
  tests::test_analytic_weight<Kernel>(spwa, idw);
  tests::test_analytic_weight<Kernel>(spwb, spwb);
}

int main(int, char**)
{
  test_kernel<SCKER>();
  test_kernel<EPICK>();
  test_kernel<EPECK>();
  std::cout << "* test_shepard_weights: SUCCESS" << std::endl;
  return EXIT_SUCCESS;
}
