#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/Delaunay_domain_2.h>
#include <CGAL/Barycentric_coordinates_2/Wachspress_coordinates_2.h>
#include <CGAL/Barycentric_coordinates_2/Discrete_harmonic_coordinates_2.h>
#include <CGAL/Barycentric_coordinates_2/Mean_value_coordinates_2.h>
#include <CGAL/Barycentric_coordinates_2/Harmonic_coordinates_2.h>

using Kernel  = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT      = typename Kernel::FT;
using Point_2 = typename Kernel::Point_2;

using Vertices = std::vector<Point_2>;
using Domain   = CGAL::Barycentric_coordinates::Delaunay_domain_2<Vertices, Kernel>;
using HMC2     = CGAL::Barycentric_coordinates::Harmonic_coordinates_2<Vertices, Domain, Kernel>;

void test_equality(
  const std::vector<FT>& coordinates_a,
  const std::vector<FT>& coordinates_b) {

  const FT epsilon = FT(1) / FT(100000);
  assert(coordinates_a.size() == 4);
  assert(coordinates_b.size() == 4);

  assert(CGAL::abs(coordinates_a[0] - coordinates_b[0]) < epsilon);
  assert(CGAL::abs(coordinates_a[1] - coordinates_b[1]) < epsilon);
  assert(CGAL::abs(coordinates_a[2] - coordinates_b[2]) < epsilon);
  assert(CGAL::abs(coordinates_a[3] - coordinates_b[3]) < epsilon);
}

int main() {

  const std::vector<Point_2> vertices = {
    Point_2(0, 0),
    Point_2(1, 0),
    Point_2(1, 1),
    Point_2(0, 1)
  };

  const FT h = FT(1) / FT(2);
  const Point_2 query(h, h);

  std::vector<FT> wp_coordinates;
  std::vector<FT> mv_coordinates;
  std::vector<FT> dh_coordinates;
  std::vector<FT> hm_coordinates;

  CGAL::Barycentric_coordinates::wachspress_coordinates_2(
    vertices, query, std::back_inserter(wp_coordinates));
  CGAL::Barycentric_coordinates::mean_value_coordinates_2(
    vertices, query, std::back_inserter(mv_coordinates));
  CGAL::Barycentric_coordinates::discrete_harmonic_coordinates_2(
    vertices, query, std::back_inserter(dh_coordinates));

  std::vector<Point_2> seeds = { query };
  Domain domain(vertices);
  domain.create(FT(1) / FT(100), seeds);
  HMC2 harmonic_coordinates_2(vertices, domain);
  harmonic_coordinates_2.compute();
  harmonic_coordinates_2(query, std::back_inserter(hm_coordinates));

  test_equality(wp_coordinates, dh_coordinates);
  test_equality(wp_coordinates, mv_coordinates);
  test_equality(wp_coordinates, hm_coordinates);
  test_equality(mv_coordinates, dh_coordinates);
  test_equality(mv_coordinates, hm_coordinates);
  test_equality(dh_coordinates, hm_coordinates);

  std::cout << "test_bc_all_coordinates: PASSED" << std::endl;
  return EXIT_SUCCESS;
}
