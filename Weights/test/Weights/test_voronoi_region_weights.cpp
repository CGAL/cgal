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

  const Point_2 p( 2, 0);
  const Point_2 q( 0, 2);
  const Point_2 r(-2, 0);
  const FT w1 = CGAL::Weights::voronoi_area(p, q, r);
  assert(w1 == FT(2));

  const Point_3 s( 0, -2, 0);
  const Point_3 t( 0,  0, 2);
  const Point_3 u( 0,  2, 0);
  const FT w4 = CGAL::Weights::voronoi_area(s, t, u);
  assert(w4 == FT(2));

  const wrappers::Voronoi_region_wrapper<Kernel> vor;
  tests::test_region_weight<Kernel>(vor);
}

int main(int, char**)
{
  test_kernel<SCKER>();
  test_kernel<EPICK>();
  test_kernel<EPECK>();
  std::cout << "* test_voronoi_region_weights: SUCCESS" << std::endl;
  return EXIT_SUCCESS;
}
