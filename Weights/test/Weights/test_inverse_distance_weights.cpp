#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include "include/utils.h"
#include "include/wrappers.h"

// Typedefs.
using SCKER = CGAL::Simple_cartesian<double>;
using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;

template<typename Kernel>
void test_overloads() {
  using FT      = typename Kernel::FT;
  using Point_2 = typename Kernel::Point_2;
  using Point_3 = typename Kernel::Point_3;
  const Point_2 p1(0, 0);
  const Point_2 q1(1, 0);
  const Point_3 p2(0, 0, 1);
  const Point_3 q2(1, 0, 1);
  const FT a2 = CGAL::Weights::inverse_distance_weight(p1, q1);
  const FT a3 = CGAL::Weights::inverse_distance_weight(p2, q2);
  assert(a2 == FT(1));
  assert(a3 == FT(1));
  assert(CGAL::Weights::inverse_distance_weight(p1, p1, q1, q1) == a2);
  assert(CGAL::Weights::inverse_distance_weight(p2, p2, q2, q2) == a3);
  struct Traits : public Kernel { };
  assert(CGAL::Weights::inverse_distance_weight(p1, p1, q1, q1, Traits()) == a2);
  assert(CGAL::Weights::inverse_distance_weight(p2, p2, q2, q2, Traits()) == a3);
}

template<typename Kernel>
bool test_kernel() {
  test_overloads<Kernel>();
  const wrappers::Inverse_distance_wrapper<Kernel> idw;
  const wrappers::Shepard_wrapper<Kernel> spw(1);
  return tests::test_analytic_weight<Kernel>(idw, spw);
}

int main() {
  assert(test_kernel<SCKER>());
  assert(test_kernel<EPICK>());
  assert(test_kernel<EPECK>());
  std::cout << "* test_inverse_distance_weights: SUCCESS" << std::endl;
  return EXIT_SUCCESS;
}
