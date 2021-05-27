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
  const FT a = CGAL::Weights::uniform_weight();
  if (a != FT(1)) return false;
  const Point_2 p(0, 0);
  const Point_3 q(0, 0, 0);
  if (CGAL::Weights::uniform_weight(p, p, p, p) != a) return false;
  if (CGAL::Weights::uniform_weight(q, q, q, q) != a) return false;
  struct Traits : public Kernel { };
  if (CGAL::Weights::uniform_weight(p, p, p, p, Traits()) != a) return false;
  if (CGAL::Weights::uniform_weight(q, q, q, q, Traits()) != a) return false;
  return true;
}

template<typename Kernel>
bool test_kernel() {
  if (!test_overloads<Kernel>()) return false;
  const wrappers::Uniform_wrapper<Kernel> uni;
  return tests::test_analytic_weight<Kernel>(uni, uni);
}

int main() {
  assert(test_kernel<SCKER>());
  assert(test_kernel<EPICK>());
  assert(test_kernel<EPECK>());
  std::cout << "* test_uniform_weights: SUCCESS" << std::endl;
  return EXIT_SUCCESS;
}
