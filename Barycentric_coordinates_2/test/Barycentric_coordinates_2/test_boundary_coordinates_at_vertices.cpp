#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/boundary_coordinates_2.h>

template<typename Kernel>
void test_boundary_coordinates_at_vertices() {

  using FT      = typename Kernel::FT;
  using Point_2 = typename Kernel::Point_2;

  const std::vector<Point_2> vertices = {
    Point_2(0, 0),
    Point_2(1, 0),
    Point_2( FT(3) / FT(2), 1),
    Point_2( FT(1) / FT(2), 2),
    Point_2(-FT(1) / FT(2), FT(3) / FT(2)),
    Point_2(-FT(1) / FT(2), FT(1) / FT(2))
  };

  const FT expected[36] = {
    1, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0,
    0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 0, 1
  };

  std::size_t count = 0;
  std::vector<FT> coordinates;
  for (std::size_t i = 0; i < 6; ++i) {
    CGAL::Barycentric_coordinates::boundary_coordinates_2(
      vertices, vertices[i], std::back_inserter(coordinates));

    assert(
      coordinates[count + 0] - expected[count + 0] == FT(0) &&
      coordinates[count + 1] - expected[count + 1] == FT(0) &&
      coordinates[count + 2] - expected[count + 2] == FT(0) &&
      coordinates[count + 3] - expected[count + 3] == FT(0) &&
      coordinates[count + 4] - expected[count + 4] == FT(0) &&
      coordinates[count + 5] - expected[count + 5] == FT(0) );
    count += 6;
  }
}

int main() {

  test_boundary_coordinates_at_vertices< CGAL::Simple_cartesian<double> >();
  test_boundary_coordinates_at_vertices< CGAL::Exact_predicates_exact_constructions_kernel >();
  test_boundary_coordinates_at_vertices< CGAL::Exact_predicates_inexact_constructions_kernel >();

  std::cout << "test_boundary_coordinates_at_vertices: PASSED" << std::endl;
  return EXIT_SUCCESS;
}
