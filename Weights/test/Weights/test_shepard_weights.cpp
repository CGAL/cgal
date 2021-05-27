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
bool test_overloads() {
  using FT      = typename Kernel::FT;
  using Point_2 = typename Kernel::Point_2;
  using Point_3 = typename Kernel::Point_3;
  const Point_2 p1(0, 0);
  const Point_2 q1(1, 0);
  const Point_3 p2(0, 0, 1);
  const Point_3 q2(1, 0, 1);
  const FT a2 = CGAL::Weights::shepard_weight(p1, q1);
  const FT a3 = CGAL::Weights::shepard_weight(p2, q2);
  if (a2 != FT(1)) return false;
  if (a3 != FT(1)) return false;
  if (CGAL::Weights::shepard_weight(p1, p1, q1, q1) != a2) return false;
  if (CGAL::Weights::shepard_weight(p2, p2, q2, q2) != a3) return false;
  struct Traits : public Kernel { };
  if (CGAL::Weights::shepard_weight(p1, p1, q1, q1, 1, Traits()) != a2) return false;
  if (CGAL::Weights::shepard_weight(p2, p2, q2, q2, 1, Traits()) != a3) return false;
  return true;
}

template<typename Kernel>
bool test_kernel() {
  if (!test_overloads<Kernel>()) return false;
  const wrappers::Shepard_wrapper<Kernel> spwa(1);
  const wrappers::Shepard_wrapper<Kernel> spwb(2);
  const wrappers::Inverse_distance_wrapper<Kernel> idw;
  if (!tests::test_analytic_weight<Kernel>(spwa, idw)) return false;
  return tests::test_analytic_weight<Kernel>(spwb, spwb);
}

int main() {
  assert(test_kernel<SCKER>());
  assert(test_kernel<EPICK>());
  assert(test_kernel<EPECK>());
  std::cout << "* test_shepard_weights: SUCCESS" << std::endl;
  return EXIT_SUCCESS;
}
