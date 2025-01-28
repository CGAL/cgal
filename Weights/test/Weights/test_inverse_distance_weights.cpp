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

  const FT w1 = CGAL::Weights::inverse_distance_weight(p1, q1);
  const FT w2 = CGAL::Weights::inverse_distance_weight(p2, q2);
  assert(w1 == FT(1) / FT(2));
  assert(w2 == FT(1) / FT(2));
  assert(CGAL::Weights::inverse_distance_weight(p1, p1, q1, q1) == w1);
  assert(CGAL::Weights::inverse_distance_weight(p2, p2, q2, q2) == w2);

  struct Traits : public Kernel { };
  assert(CGAL::Weights::inverse_distance_weight(p1, q1, Traits()) == w1);
  assert(CGAL::Weights::inverse_distance_weight(p2, q2, Traits()) == w2);
  assert(CGAL::Weights::inverse_distance_weight(p1, p1, q1, q1, Traits()) == w1);
  assert(CGAL::Weights::inverse_distance_weight(p2, p2, q2, q2, Traits()) == w2);
}

template<typename Kernel>
void test_kernel()
{
  test_overloads<Kernel>();
  const wrappers::Inverse_distance_wrapper<Kernel> idw;
  const wrappers::Shepard_wrapper<Kernel> spw(1);
  tests::test_analytic_weight<Kernel>(idw, spw);
}

int main(int, char**)
{
  test_kernel<SCKER>();
  test_kernel<EPICK>();
  test_kernel<EPECK>();
  std::cout << "* test_inverse_distance_weights: SUCCESS" << std::endl;
  return EXIT_SUCCESS;
}
