#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/boundary_coordinates_2.h>

template<typename Kernel>
void test_boundary_coordinates_on_edges() {

  using FT      = typename Kernel::FT;
  using Point_2 = typename Kernel::Point_2;

  const std::vector<Point_2> vertices = {
    Point_2(0, 0),
    Point_2(1, 0),
    Point_2( FT(1) / FT(2), 1),
    Point_2( FT(3) / FT(2), FT(3) / FT(2)),
    Point_2(-FT(1) / FT(2), FT(3) / FT(2)),
    Point_2(0, 1),
    Point_2(-FT(1) / FT(2), FT(1) / FT(2))
  };

  const Point_2 queries[7] = {
    Point_2(FT(1) / FT(2), 0),
    Point_2(FT(3) / FT(4), FT(1) / FT(2)),
    Point_2(FT(1) / FT(2), 1),
    Point_2(FT(3) / FT(4), FT(3) / FT(2)),
    Point_2(0, 1),
    Point_2(-FT(1) / FT(8), FT(7) / FT(8)),
    Point_2(-FT(1) / FT(8), FT(1) / FT(8))
  };

  const FT expected[49] = {
    FT(1) / FT(2), FT(1) / FT(2), 0, 0, 0, 0, 0,
    0, FT(1) / FT(2), FT(1) / FT(2), 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, FT(5) / FT(8), FT(3) / FT(8), 0, 0,
    0, 0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 0, FT(3) / FT(4), FT(1) / FT(4),
    FT(3) / FT(4), 0, 0, 0, 0, 0, FT(1) / FT(4)
  };

  std::size_t count = 0;
  std::vector<FT> coordinates;

  CGAL::Barycentric_coordinates::boundary_coordinates_2(
    vertices, vertices[2], std::back_inserter(coordinates));
  assert(
    coordinates[count + 0] - expected[14 + 0] == FT(0) &&
    coordinates[count + 1] - expected[14 + 1] == FT(0) &&
    coordinates[count + 2] - expected[14 + 2] == FT(0) &&
    coordinates[count + 3] - expected[14 + 3] == FT(0) &&
    coordinates[count + 4] - expected[14 + 4] == FT(0) &&
    coordinates[count + 5] - expected[14 + 5] == FT(0) &&
    coordinates[count + 6] - expected[14 + 6] == FT(0) );
  count += 7;

  CGAL::Barycentric_coordinates::boundary_coordinates_2(
    vertices, vertices[5], std::back_inserter(coordinates));
  assert(
    coordinates[count + 0] - expected[28 + 0] == FT(0) &&
    coordinates[count + 1] - expected[28 + 1] == FT(0) &&
    coordinates[count + 2] - expected[28 + 2] == FT(0) &&
    coordinates[count + 3] - expected[28 + 3] == FT(0) &&
    coordinates[count + 4] - expected[28 + 4] == FT(0) &&
    coordinates[count + 5] - expected[28 + 5] == FT(0) &&
    coordinates[count + 6] - expected[28 + 6] == FT(0) );
  count += 7;

  CGAL::Barycentric_coordinates::boundary_coordinates_2(
    vertices, queries[2], std::back_inserter(coordinates));
  assert(
    coordinates[count + 0] - expected[14 + 0] == FT(0) &&
    coordinates[count + 1] - expected[14 + 1] == FT(0) &&
    coordinates[count + 2] - expected[14 + 2] == FT(0) &&
    coordinates[count + 3] - expected[14 + 3] == FT(0) &&
    coordinates[count + 4] - expected[14 + 4] == FT(0) &&
    coordinates[count + 5] - expected[14 + 5] == FT(0) &&
    coordinates[count + 6] - expected[14 + 6] == FT(0) );
  count += 7;

  CGAL::Barycentric_coordinates::boundary_coordinates_2(
    vertices, queries[4], std::back_inserter(coordinates));
  assert(
    coordinates[count + 0] - expected[28 + 0] == FT(0) &&
    coordinates[count + 1] - expected[28 + 1] == FT(0) &&
    coordinates[count + 2] - expected[28 + 2] == FT(0) &&
    coordinates[count + 3] - expected[28 + 3] == FT(0) &&
    coordinates[count + 4] - expected[28 + 4] == FT(0) &&
    coordinates[count + 5] - expected[28 + 5] == FT(0) &&
    coordinates[count + 6] - expected[28 + 6] == FT(0) );
  count += 7;

  count = 0;
  coordinates.clear();
  for (std::size_t i = 0; i < 7; ++i) {

    CGAL::Barycentric_coordinates::boundary_coordinates_2(
      vertices, queries[i], std::back_inserter(coordinates));
    assert(
      coordinates[count + 0] - expected[count + 0] == FT(0) &&
      coordinates[count + 1] - expected[count + 1] == FT(0) &&
      coordinates[count + 2] - expected[count + 2] == FT(0) &&
      coordinates[count + 3] - expected[count + 3] == FT(0) &&
      coordinates[count + 4] - expected[count + 4] == FT(0) &&
      coordinates[count + 5] - expected[count + 5] == FT(0) &&
      coordinates[count + 6] - expected[count + 6] == FT(0) );
    count += 7;
  }
}

int main() {

  test_boundary_coordinates_on_edges< CGAL::Simple_cartesian<double> >();
  test_boundary_coordinates_on_edges< CGAL::Exact_predicates_exact_constructions_kernel >();
  test_boundary_coordinates_on_edges< CGAL::Exact_predicates_inexact_constructions_kernel >();

  std::cout << "test_boundary_coordinates_on_edges: PASSED" << std::endl;
  return EXIT_SUCCESS;
}
