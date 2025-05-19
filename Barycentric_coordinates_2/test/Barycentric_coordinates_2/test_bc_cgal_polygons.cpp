#include <CGAL/Polygon_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/Wachspress_coordinates_2.h>
#include <CGAL/Barycentric_coordinates_2/Discrete_harmonic_coordinates_2.h>
#include <CGAL/Barycentric_coordinates_2/Mean_value_coordinates_2.h>
#include <CGAL/Barycentric_coordinates_2/boundary_coordinates_2.h>

using Kernel    = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT        = typename Kernel::FT;
using Point_2   = typename Kernel::Point_2;
using Polygon_2 = CGAL::Polygon_2<Kernel>;

void check_result(
  const std::vector<FT>& weights,
  const std::vector<FT>& coordinates) {

  assert(coordinates.size() == 4);
  assert(weights.size() == coordinates.size());

  const FT quater = FT(1) / FT(4);
  const FT epsilon = FT(1) / FT(100000000000000);
  assert(
    CGAL::abs(coordinates[0] - quater) < epsilon &&
    CGAL::abs(coordinates[1] - quater) < epsilon &&
    CGAL::abs(coordinates[2] - quater) < epsilon &&
    CGAL::abs(coordinates[3] - quater) < epsilon );

  const FT sum = weights[0] + weights[1] + weights[2] + weights[3];
  assert(
    CGAL::abs(coordinates[0] - weights[0] / sum) < epsilon &&
    CGAL::abs(coordinates[1] - weights[1] / sum) < epsilon &&
    CGAL::abs(coordinates[2] - weights[2] / sum) < epsilon &&
    CGAL::abs(coordinates[3] - weights[3] / sum) < epsilon );
}

int main() {

  Polygon_2 polygon;
  polygon.push_back(Point_2(0, 0));
  polygon.push_back(Point_2(1, 0));
  polygon.push_back(Point_2(1, 1));
  polygon.push_back(Point_2(0, 1));

  const FT h = FT(1) / FT(2);
  const Point_2 query(h, h);

  std::vector<FT> weights;
  weights.reserve(4);
  std::vector<FT> coordinates;
  coordinates.reserve(4);

  weights.clear(); coordinates.clear();
  CGAL::Barycentric_coordinates::wachspress_weights_2(
    polygon, query, std::back_inserter(weights));
  CGAL::Barycentric_coordinates::wachspress_coordinates_2(
    polygon, query, std::back_inserter(coordinates));
  check_result(weights, coordinates);

  weights.clear(); coordinates.clear();
  CGAL::Barycentric_coordinates::mean_value_weights_2(
    polygon, query, std::back_inserter(weights));
  CGAL::Barycentric_coordinates::mean_value_coordinates_2(
    polygon, query, std::back_inserter(coordinates));
  check_result(weights, coordinates);

  weights.clear(); coordinates.clear();
  CGAL::Barycentric_coordinates::discrete_harmonic_weights_2(
    polygon, query, std::back_inserter(weights));
  CGAL::Barycentric_coordinates::discrete_harmonic_coordinates_2(
    polygon, query, std::back_inserter(coordinates));
  check_result(weights, coordinates);

  coordinates.clear();
  CGAL::Barycentric_coordinates::boundary_coordinates_2(
    polygon, query, std::back_inserter(coordinates));
  assert(
    coordinates[0] == FT(0) && coordinates[1] == FT(0) &&
    coordinates[2] == FT(0) && coordinates[3] == FT(0) );

  std::cout << "test_bc_cgal_polygons: PASSED" << std::endl;
  return EXIT_SUCCESS;
}
