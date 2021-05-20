#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/triangle_coordinates_2.h>

template<typename Kernel>
void test_triangle_coordinates() {

  using FT      = typename Kernel::FT;
  using Point_2 = typename Kernel::Point_2;

  const Point_2 p0 = Point_2(0, 0);
  const Point_2 p1 = Point_2(1, 0);
  const Point_2 p2 = Point_2(FT(1) / FT(2), 1);

  const Point_2 queries[5] = {
    Point_2(FT(1) / FT(2), FT(1) / FT(2)),
    Point_2(FT(1) / FT(2), 0),
    Point_2(FT(1) / FT(2), 1),
    Point_2(1, FT(1) / FT(2)),
    Point_2(0, FT(1) / FT(2))
  };

  const FT expected[15] = {
    FT(1) / FT(4), FT(1) / FT(4), FT(1) / FT(2),
    FT(1) / FT(2), FT(1) / FT(2), 0,
    0, 0, 1,
   -FT(1) / FT(4),  FT(3) / FT(4), FT(1) / FT(2),
    FT(3) / FT(4), -FT(1) / FT(4), FT(1) / FT(2)
  };

  std::size_t count = 0;
  std::vector<FT> coordinates;
  for (std::size_t i = 0; i < 5; ++i) {
    CGAL::Barycentric_coordinates::triangle_coordinates_2(
      p0, p1, p2, queries[i], std::back_inserter(coordinates));
    assert(
      coordinates[count + 0] - expected[count + 0] == FT(0) &&
      coordinates[count + 1] - expected[count + 1] == FT(0) &&
      coordinates[count + 2] - expected[count + 2] == FT(0) );
    count += 3;
  }

  count = 0;
  for (std::size_t i = 0; i < 5; ++i) {
    const auto tuple =
      CGAL::Barycentric_coordinates::triangle_coordinates_in_tuple_2(
        p0, p1, p2, queries[i]);
    assert(
      std::get<0>(tuple) - expected[count + 0] == FT(0) &&
      std::get<1>(tuple) - expected[count + 1] == FT(0) &&
      std::get<2>(tuple) - expected[count + 2] == FT(0) );
    count += 3;
  }
}

int main() {

  test_triangle_coordinates< CGAL::Simple_cartesian<double> >();
  test_triangle_coordinates< CGAL::Exact_predicates_exact_constructions_kernel >();
  test_triangle_coordinates< CGAL::Exact_predicates_inexact_constructions_kernel >();

  std::cout << "test_triangle_coordinates: PASSED" << std::endl;
  return EXIT_SUCCESS;
}
