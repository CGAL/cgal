#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Projection_traits_xz_3.h>
#include <CGAL/Projection_traits_yz_3.h>
#include <CGAL/Weights.h>

// Typedefs.
using SCKER = CGAL::Simple_cartesian<double>;
using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;

template<typename Kernel>
void test_kernel() {
  using FT      = typename Kernel::FT;
  using Point_2 = typename Kernel::Point_2;
  using Point_3 = typename Kernel::Point_3;

  using XY_Traits = CGAL::Projection_traits_xy_3<Kernel>;
  using XZ_Traits = CGAL::Projection_traits_xz_3<Kernel>;
  using YZ_Traits = CGAL::Projection_traits_yz_3<Kernel>;

  const XY_Traits xy_traits;
  const XZ_Traits xz_traits;
  const YZ_Traits yz_traits;

  const Point_2 t(-1,  0);
  const Point_2 r( 0, -1);
  const Point_2 p( 1,  0);
  const Point_2 q( 0,  0);

  // XY.
  const Point_3 t1(t.x(), t.y(), 1);
  const Point_3 r1(r.x(), r.y(), 1);
  const Point_3 p1(p.x(), p.y(), 1);
  const Point_3 q1(q.x(), q.y(), 1);

  // XZ.
  const Point_3 t2(t.x(), 1, t.y());
  const Point_3 r2(r.x(), 1, r.y());
  const Point_3 p2(p.x(), 1, p.y());
  const Point_3 q2(q.x(), 1, q.y());

  // YZ.
  const Point_3 t3(1, t.x(), t.y());
  const Point_3 r3(1, r.x(), r.y());
  const Point_3 p3(1, p.x(), p.y());
  const Point_3 q3(1, q.x(), q.y());

  const FT ref_value = FT(4);
  assert(CGAL::Weights::authalic_weight(t1, r1, p1, q1, xy_traits) == ref_value);
  assert(CGAL::Weights::authalic_weight(t2, r2, p2, q2, xz_traits) == ref_value);
  assert(CGAL::Weights::authalic_weight(t3, r3, p3, q3, yz_traits) == ref_value);

  assert(CGAL::Weights::wachspress_weight(t1, r1, p1, q1, xy_traits) == ref_value);
  assert(CGAL::Weights::wachspress_weight(t2, r2, p2, q2, xz_traits) == ref_value);
  assert(CGAL::Weights::wachspress_weight(t3, r3, p3, q3, yz_traits) == ref_value);

  assert(CGAL::Weights::cotangent_weight(t1, r1, p1, q1, xy_traits) == ref_value);
  assert(CGAL::Weights::cotangent_weight(t2, r2, p2, q2, xz_traits) == ref_value);
  assert(CGAL::Weights::cotangent_weight(t3, r3, p3, q3, yz_traits) == ref_value);

  assert(CGAL::Weights::discrete_harmonic_weight(t1, r1, p1, q1, xy_traits) == ref_value);
  assert(CGAL::Weights::discrete_harmonic_weight(t2, r2, p2, q2, xz_traits) == ref_value);
  assert(CGAL::Weights::discrete_harmonic_weight(t3, r3, p3, q3, yz_traits) == ref_value);

  assert(CGAL::Weights::tangent_weight(t1, r1, p1, q1, xy_traits) == ref_value);
  assert(CGAL::Weights::tangent_weight(t2, r2, p2, q2, xz_traits) == ref_value);
  assert(CGAL::Weights::tangent_weight(t3, r3, p3, q3, yz_traits) == ref_value);

  assert(CGAL::Weights::mean_value_weight(t1, r1, p1, q1, xy_traits) == ref_value);
  assert(CGAL::Weights::mean_value_weight(t2, r2, p2, q2, xz_traits) == ref_value);
  assert(CGAL::Weights::mean_value_weight(t3, r3, p3, q3, yz_traits) == ref_value);

  assert(CGAL::Weights::uniform_area(t1, r1, p1, xy_traits) == FT(1));
  assert(CGAL::Weights::uniform_area(t2, r2, p2, xz_traits) == FT(1));
  assert(CGAL::Weights::uniform_area(t3, r3, p3, yz_traits) == FT(1));

  assert(CGAL::Weights::triangular_area(t1, r1, p1, xy_traits) >= FT(0));
  assert(CGAL::Weights::triangular_area(t2, r2, p2, xz_traits) >= FT(0));
  assert(CGAL::Weights::triangular_area(t3, r3, p3, yz_traits) >= FT(0));

  assert(CGAL::Weights::barycentric_area(t1, r1, p1, xy_traits) >= FT(0));
  assert(CGAL::Weights::barycentric_area(t2, r2, p2, xz_traits) >= FT(0));
  assert(CGAL::Weights::barycentric_area(t3, r3, p3, yz_traits) >= FT(0));

  assert(CGAL::Weights::voronoi_area(t1, r1, p1, xy_traits) >= FT(0));
  assert(CGAL::Weights::voronoi_area(t2, r2, p2, xz_traits) >= FT(0));
  assert(CGAL::Weights::voronoi_area(t3, r3, p3, yz_traits) >= FT(0));

  assert(CGAL::Weights::mixed_voronoi_area(t1, r1, p1, xy_traits) >= FT(0));
  assert(CGAL::Weights::mixed_voronoi_area(t2, r2, p2, xz_traits) >= FT(0));
  assert(CGAL::Weights::mixed_voronoi_area(t3, r3, p3, yz_traits) >= FT(0));

  assert(CGAL::Weights::uniform_weight(t1, r1, p1, q1, xy_traits) == FT(1));
  assert(CGAL::Weights::uniform_weight(t2, r2, p2, q2, xz_traits) == FT(1));
  assert(CGAL::Weights::uniform_weight(t3, r3, p3, q3, yz_traits) == FT(1));

  assert(CGAL::Weights::inverse_distance_weight(t1, r1, p1, q1, xy_traits) == FT(1));
  assert(CGAL::Weights::inverse_distance_weight(t2, r2, p2, q2, xz_traits) == FT(1));
  assert(CGAL::Weights::inverse_distance_weight(t3, r3, p3, q3, yz_traits) == FT(1));

  assert(CGAL::Weights::shepard_weight(t1, r1, p1, q1, 1, xy_traits) == FT(1));
  assert(CGAL::Weights::shepard_weight(t2, r2, p2, q2, 1, xz_traits) == FT(1));
  assert(CGAL::Weights::shepard_weight(t3, r3, p3, q3, 1, yz_traits) == FT(1));

  assert(CGAL::Weights::three_point_family_weight(t1, r1, p1, q1, 1, xy_traits) == ref_value);
  assert(CGAL::Weights::three_point_family_weight(t2, r2, p2, q2, 1, xz_traits) == ref_value);
  assert(CGAL::Weights::three_point_family_weight(t3, r3, p3, q3, 1, yz_traits) == ref_value);
}

int main() {
  test_kernel<SCKER>();
  test_kernel<EPICK>();
  test_kernel<EPECK>();
  std::cout << "* test_projected_weights: SUCCESS" << std::endl;
  return EXIT_SUCCESS;
}
