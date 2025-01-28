#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Projection_traits_xz_3.h>
#include <CGAL/Projection_traits_yz_3.h>
#include <CGAL/Barycentric_coordinates_2/Delaunay_domain_2.h>
#include <CGAL/Barycentric_coordinates_2/Wachspress_coordinates_2.h>
#include <CGAL/Barycentric_coordinates_2/Discrete_harmonic_coordinates_2.h>
#include <CGAL/Barycentric_coordinates_2/Mean_value_coordinates_2.h>
#include <CGAL/Barycentric_coordinates_2/Harmonic_coordinates_2.h>
#include <CGAL/Barycentric_coordinates_2/boundary_coordinates_2.h>

using Kernel  = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT      = typename Kernel::FT;
using Point_3 = typename Kernel::Point_3;

using Projection_traits_xy = CGAL::Projection_traits_xy_3<Kernel>;
using Projection_traits_xz = CGAL::Projection_traits_xz_3<Kernel>;
using Projection_traits_yz = CGAL::Projection_traits_yz_3<Kernel>;

using Vertices  = std::vector<Point_3>;
using Point_map = CGAL::Identity_property_map<Point_3>;

void check_result(
  const std::vector<FT>& weights,
  const std::vector<FT>& coordinates,
  const FT epsilon = FT(1) / FT(100000000000000)) {

  const FT quater = FT(1) / FT(4);
  assert(coordinates.size() == 4);
  assert(
    CGAL::abs(coordinates[0] - quater) < epsilon &&
    CGAL::abs(coordinates[1] - quater) < epsilon &&
    CGAL::abs(coordinates[2] - quater) < epsilon &&
    CGAL::abs(coordinates[3] - quater) < epsilon );

  if (weights.size() == 4) {
    assert(coordinates.size() == weights.size());
    const FT sum = weights[0] + weights[1] + weights[2] + weights[3];
    assert(
      CGAL::abs(coordinates[0] - weights[0] / sum) < epsilon &&
      CGAL::abs(coordinates[1] - weights[1] / sum) < epsilon &&
      CGAL::abs(coordinates[2] - weights[2] / sum) < epsilon &&
      CGAL::abs(coordinates[3] - weights[3] / sum) < epsilon );
  }
}

template<typename Projection_traits>
void test_projection_traits(
  const std::vector<Point_3>& vertices,
  const Point_3& query,
  const Projection_traits& projection_traits) {

  std::vector<FT> weights;
  std::vector<FT> coordinates;

  weights.clear(); coordinates.clear();
  CGAL::Barycentric_coordinates::wachspress_weights_2(
    vertices, query, std::back_inserter(weights), projection_traits);
  CGAL::Barycentric_coordinates::wachspress_coordinates_2(
    vertices, query, std::back_inserter(coordinates), projection_traits);
  check_result(weights, coordinates);

  weights.clear(); coordinates.clear();
  CGAL::Barycentric_coordinates::mean_value_weights_2(
    vertices, query, std::back_inserter(weights), projection_traits);
  CGAL::Barycentric_coordinates::mean_value_coordinates_2(
    vertices, query, std::back_inserter(coordinates), projection_traits);
  check_result(weights, coordinates);

  weights.clear(); coordinates.clear();
  CGAL::Barycentric_coordinates::discrete_harmonic_weights_2(
    vertices, query, std::back_inserter(weights), projection_traits);
  CGAL::Barycentric_coordinates::discrete_harmonic_coordinates_2(
    vertices, query, std::back_inserter(coordinates), projection_traits);
  check_result(weights, coordinates);

  using Domain = CGAL::Barycentric_coordinates::
    Delaunay_domain_2<Vertices, Projection_traits, Point_map>;
  using HMC2   = CGAL::Barycentric_coordinates::
    Harmonic_coordinates_2<Vertices, Domain, Projection_traits, Point_map>;

  coordinates.clear();
  Point_map point_map;

  CGAL::Barycentric_coordinates::boundary_coordinates_2(
    vertices, vertices[0], std::back_inserter(coordinates), projection_traits, point_map);
  assert(
    coordinates[0] == FT(1) && coordinates[1] == FT(0) &&
    coordinates[2] == FT(0) && coordinates[3] == FT(0) );

  const FT max_edge_length = FT(1) / FT(10);
  const std::vector<Point_3> seeds = { query };
  Domain domain(vertices, projection_traits, point_map);
  domain.create(max_edge_length, seeds);

  HMC2 harmonic_coordinates_2(
    vertices, domain, projection_traits, point_map);
  harmonic_coordinates_2.compute();
  coordinates.clear();
  for (std::size_t k = 0; k < domain.number_of_vertices(); ++k)
    harmonic_coordinates_2(k, std::back_inserter(coordinates));
  assert(coordinates.size() ==
    domain.number_of_vertices() * vertices.size());

  weights.clear(); coordinates.clear();
  harmonic_coordinates_2(query, std::back_inserter(coordinates));
  check_result(weights, coordinates, max_edge_length);
}

int main() {
  const FT h = FT(1) / FT(2);

  // Test XY plane.
  Projection_traits_xy projection_traits_xy;
  const std::vector<Point_3> vertices_xy = {
    Point_3(0, 0, 1),
    Point_3(1, 0, 1),
    Point_3(1, 1, 1),
    Point_3(0, 1, 1)
  };
  const Point_3 query_xy(h, h, 1);
  test_projection_traits(vertices_xy, query_xy, projection_traits_xy);

  // Test XZ plane.
  Projection_traits_xz projection_traits_xz;
  const std::vector<Point_3> vertices_xz = {
    Point_3(0, 1, 0),
    Point_3(1, 1, 0),
    Point_3(1, 1, 1),
    Point_3(0, 1, 1)
  };
  const Point_3 query_xz(h, 1, h);
  test_projection_traits(vertices_xz, query_xz, projection_traits_xz);

  // Test YZ plane.
  Projection_traits_yz projection_traits_yz;
  const std::vector<Point_3> vertices_yz = {
    Point_3(1, 0, 0),
    Point_3(1, 1, 0),
    Point_3(1, 1, 1),
    Point_3(1, 0, 1)
  };
  const Point_3 query_yz(1, h, h);
  test_projection_traits(vertices_yz, query_yz, projection_traits_yz);

  std::cout << "test_bc_projection_traits: PASSED" << std::endl;
  return EXIT_SUCCESS;
}
