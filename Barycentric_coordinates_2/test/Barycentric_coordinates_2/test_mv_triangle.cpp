#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/triangle_coordinates_2.h>
#include <CGAL/Barycentric_coordinates_2/Mean_value_coordinates_2.h>

template<typename Kernel>
void test_mv_triangle() {

  using FT      = typename Kernel::FT;
  using Point_2 = typename Kernel::Point_2;

  const std::vector<Point_2> vertices = {
    Point_2(0, 0),
    Point_2(1, 0),
    Point_2(0, 1)
  };

  std::vector<FT> tri_coordinates;
  std::vector<FT>  mv_coordinates;

  const FT step  = FT(1) / FT(100);
  const FT scale = FT(50);

  std::size_t count = 0;
  const FT limit = scale * step;
  const FT epsilon = FT(1) / FT(100000000000000);

  for (FT x = step; x < limit; x += step) {
    for (FT y = step; y < limit; y += step) {
      const Point_2 query(x, y);

      CGAL::Barycentric_coordinates::triangle_coordinates_2(
        vertices[0], vertices[1], vertices[2], query, std::back_inserter(tri_coordinates));
      CGAL::Barycentric_coordinates::mean_value_coordinates_2(
        vertices, query, std::back_inserter(mv_coordinates));

      assert(tri_coordinates[count + 0] >= FT(0) && tri_coordinates[count + 0] <= FT(1));
      assert(tri_coordinates[count + 1] >= FT(0) && tri_coordinates[count + 1] <= FT(1));
      assert(tri_coordinates[count + 2] >= FT(0) && tri_coordinates[count + 2] <= FT(1));

      assert(mv_coordinates[count + 0] >= FT(0) && mv_coordinates[count + 0] <= FT(1));
      assert(mv_coordinates[count + 1] >= FT(0) && mv_coordinates[count + 1] <= FT(1));
      assert(mv_coordinates[count + 2] >= FT(0) && mv_coordinates[count + 2] <= FT(1));

      assert(
        CGAL::abs(tri_coordinates[count + 0] - mv_coordinates[count + 0]) < epsilon &&
        CGAL::abs(tri_coordinates[count + 1] - mv_coordinates[count + 1]) < epsilon &&
        CGAL::abs(tri_coordinates[count + 2] - mv_coordinates[count + 2]) < epsilon );
      count += 3;
    }
  }
}

int main() {

  test_mv_triangle< CGAL::Simple_cartesian<double> >();
  test_mv_triangle< CGAL::Exact_predicates_inexact_constructions_kernel >();

  std::cout << "test_mv_triangle: PASSED" << std::endl;
  return EXIT_SUCCESS;
}
