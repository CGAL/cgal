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
  const Point_2 t1(-1,  0);
  const Point_2 r1( 0, -1);
  const Point_2 p1( 1,  0);
  const Point_2 q1( 0,  0);
  const Point_3 t2(-1,  0, 1);
  const Point_3 r2( 0, -1, 1);
  const Point_3 p2( 1,  0, 1);
  const Point_3 q2( 0,  0, 1);
  const FT a2 = CGAL::Weights::mean_value_weight(t1, r1, p1, q1);
  const FT a3 = CGAL::Weights::internal::mean_value_weight(t2, r2, p2, q2);
  assert(a2 >= FT(0));
  assert(a3 >= FT(0));
  assert(a2 == a3);
  struct Traits : public Kernel { };
  assert(CGAL::Weights::mean_value_weight(t1, r1, p1, q1, Traits()) == a2);
  assert(CGAL::Weights::internal::mean_value_weight(t2, r2, p2, q2, Traits()) == a3);
  CGAL::Projection_traits_xy_3<Kernel> ptraits;
  const FT a23 = CGAL::Weights::mean_value_weight(t2, r2, p2, q2, ptraits);
  assert(a23 >= FT(0));
  assert(a23 == a2 && a23 == a3);
}

template<typename Kernel>
bool test_kernel() {
  test_overloads<Kernel>();
  const wrappers::Mean_value_wrapper<Kernel> mvw;
  const wrappers::Tangent_wrapper<Kernel> tan;
  return tests::test_barycentric_weight<Kernel>(mvw, tan);
}

int main() {
  assert(test_kernel<SCKER>());
  assert(test_kernel<EPICK>());
  assert(test_kernel<EPECK>());
  std::cout << "* test_mean_value_weights: SUCCESS" << std::endl;
  return EXIT_SUCCESS;
}
