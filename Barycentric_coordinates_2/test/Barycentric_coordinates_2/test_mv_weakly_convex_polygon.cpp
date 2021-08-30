#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/Mean_value_coordinates_2.h>

template<typename Kernel>
void test_mv_weakly_convex_polygon() {

  using FT      = typename Kernel::FT;
  using Point_2 = typename Kernel::Point_2;

  const std::vector<Point_2> vertices = {
    Point_2(0, 0),
    Point_2(1, 0),
    Point_2(1, FT(1) / FT(2)),
    Point_2(1, 1),
    Point_2(0, 1),
    Point_2(0, FT(1) / FT(2)),
  };

  const Point_2 queries[7] = {
    Point_2(FT(1) - (FT(1) / FT(std::pow(10.0, 300.0))), FT(1) / FT(4)),
    Point_2(FT(1) - (FT(1) / FT(std::pow(10.0, 300.0))), FT(5) / FT(8)),
    Point_2(FT(1) - (FT(1) / FT(std::pow(10.0, 300.0))), FT(7) / FT(8)),
    Point_2((FT(1) / FT(std::pow(10.0, 300.0))), FT(5) / FT(8)),
    Point_2((FT(1) / FT(std::pow(10.0, 300.0))), FT(3) / FT(8)),
    Point_2((FT(1) / FT(std::pow(10.0, 300.0))), FT(1) / FT(4)),
    Point_2((FT(1) / FT(std::pow(10.0, 300.0))), (FT(1) / FT(std::pow(10.0, 300.0))))
  };

  const FT expected[42] = {
    0, FT(1) / FT(2), FT(1) / FT(2), 0, 0, 0,
    0, 0, FT(3) / FT(4), FT(1) / FT(4), 0, 0,
    0, 0, FT(1) / FT(4), FT(3) / FT(4), 0, 0,
    0, 0, 0, 0, FT(1) / FT(4), FT(3) / FT(4),
    FT(1) / FT(4), 0, 0, 0, 0, FT(3) / FT(4),
    FT(3) / FT(4), 0, 0, 0, 0, FT(1) / FT(4),
    1, 0, 0, 0, 0, 0
  };

  std::size_t count = 0;
  const FT epsilon = FT(1) / FT(1000000000000000);

  std::vector<FT> coordinates;
  for (std::size_t i = 0; i < 5; ++i) {
    CGAL::Barycentric_coordinates::mean_value_coordinates_2(
      vertices, queries[i], std::back_inserter(coordinates));
    assert(
      CGAL::abs(coordinates[count + 0] - expected[count + 0]) < epsilon &&
      CGAL::abs(coordinates[count + 1] - expected[count + 1]) < epsilon &&
      CGAL::abs(coordinates[count + 2] - expected[count + 2]) < epsilon &&
      CGAL::abs(coordinates[count + 3] - expected[count + 3]) < epsilon &&
      CGAL::abs(coordinates[count + 4] - expected[count + 4]) < epsilon &&
      CGAL::abs(coordinates[count + 5] - expected[count + 5]) < epsilon );
    count += 6;
  }
}

int main() {

  test_mv_weakly_convex_polygon< CGAL::Simple_cartesian<double> >();
  test_mv_weakly_convex_polygon< CGAL::Exact_predicates_inexact_constructions_kernel >();

  std::cout << "test_mv_weakly_convex_polygon: PASSED" << std::endl;
  return EXIT_SUCCESS;
}
