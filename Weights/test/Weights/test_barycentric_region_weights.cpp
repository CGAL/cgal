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
  using Point_3 = typename Kernel::Point_3;

  const Point_2 p( 1, 0);
  const Point_2 q( 0, 6);
  const Point_2 r(-1, 0);
  const FT w1 = CGAL::Weights::barycentric_area(p, q, r);
  const FT w2 = CGAL::Weights::barycentric_area(r, p, q);
  const FT w3 = CGAL::Weights::barycentric_area(q, r, p);
  assert(w1 == FT(2)); // medians subdivide a triangle into 6 triangles of equal areas
  assert(w1 == w2 && w2 == w3);

  const Point_3 s( 0, -1, 0);
  const Point_3 t( 0,  0, 6);
  const Point_3 u( 0,  1, 0);
  const FT w4 = CGAL::Weights::barycentric_area(s, t, u);
  const FT w5 = CGAL::Weights::barycentric_area(t, u, s);
  const FT w6 = CGAL::Weights::barycentric_area(u, s, t);
  assert(w4 == FT(2));
  assert(w4 == w5 && w5 == w6);

  const wrappers::Barycentric_region_wrapper<Kernel> bar;
  tests::test_region_weight<Kernel>(bar);
}

int main(int, char**)
{
  test_kernel<SCKER>();
  test_kernel<EPICK>();
  test_kernel<EPECK>();
  std::cout << "* test_barycentric_region_weights: SUCCESS" << std::endl;
  return EXIT_SUCCESS;
}
