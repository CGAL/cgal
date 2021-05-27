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
  const Point_2 t1(-1,  0);
  const Point_2 r1( 0, -1);
  const Point_2 p1( 1,  0);
  const Point_2 q1( 0,  0);
  const Point_3 t2(-1,  0, 1);
  const Point_3 r2( 0, -1, 1);
  const Point_3 p2( 1,  0, 1);
  const Point_3 q2( 0,  0, 1);
  const FT a2 = CGAL::Weights::three_point_family_weight(t1, r1, p1, q1);
  const FT a3 = CGAL::Weights::internal::three_point_family_weight(t2, r2, p2, q2);
  if (a2 < FT(0)) return false;
  if (a3 < FT(0)) return false;
  if (a2 != a3)   return false;
  struct Traits : public Kernel { };
  if (CGAL::Weights::three_point_family_weight(t1, r1, p1, q1, 1, Traits()) != a2) return false;
  if (CGAL::Weights::internal::three_point_family_weight(t2, r2, p2, q2, 1, Traits()) != a3) return false;
  return true;
}

template<typename Kernel>
bool test_kernel() {
  if (!test_overloads<Kernel>()) return false;
  using FT = typename Kernel::FT;
  const FT h = FT(1) / FT(2);
  const wrappers::Three_point_family_wrapper<Kernel> tpfa(0);
  const wrappers::Three_point_family_wrapper<Kernel> tpfb(1);
  const wrappers::Three_point_family_wrapper<Kernel> tpfc(2);
  const wrappers::Three_point_family_wrapper<Kernel> tpfd(h);
  const wrappers::Wachspress_wrapper<Kernel> whp;
  const wrappers::Mean_value_wrapper<Kernel> mvw;
  const wrappers::Discrete_harmonic_wrapper<Kernel> dhw;
  if (!tests::test_analytic_weight<Kernel>(tpfa, whp)) return false;
  if (!tests::test_analytic_weight<Kernel>(tpfb, mvw)) return false;
  if (!tests::test_analytic_weight<Kernel>(tpfc, dhw)) return false;
  return tests::test_analytic_weight<Kernel>(tpfd, tpfd);
}

int main() {
  assert(test_kernel<SCKER>());
  assert(test_kernel<EPICK>());
  assert(test_kernel<EPECK>());
  std::cout << "* test_three_point_family_weights: SUCCESS" << std::endl;
  return EXIT_SUCCESS;
}
