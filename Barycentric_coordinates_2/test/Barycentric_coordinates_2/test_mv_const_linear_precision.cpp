#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/Mean_value_coordinates_2.h>

template<typename Kernel>
void test_mv_const_linear_precision() {

  using FT      = typename Kernel::FT;
  using Point_2 = typename Kernel::Point_2;

  const std::vector<Point_2> vertices = {
    Point_2(0, 0),
    Point_2(FT(3)  / FT(4),  FT(1) / FT(2)),
    Point_2(FT(7)  / FT(4), -FT(1) / FT(2)),
    Point_2(FT(11) / FT(4),  FT(1) / FT(2)),
    Point_2(FT(7)  / FT(2), 0),
    Point_2(FT(7)  / FT(2), 2),
    Point_2(FT(11) / FT(4), FT(3) / FT(2)),
    Point_2(FT(7)  / FT(4), FT(5) / FT(2)),
    Point_2(FT(3)  / FT(4), FT(3) / FT(2)),
    Point_2(0, 2)
  };

  const FT step = FT(1) / FT(100);
  const FT x_scale = FT(350);
  const FT y_scale = FT(100);

  std::size_t count = 0;
  const FT epsilon = FT(1) / FT(100000000);

  const FT limit_x = x_scale * step;
  const FT half    = FT(1) / FT(2);
  const FT y_start = half + step;
  const FT limit_y = half + y_scale * step;

  std::vector<FT> coordinates;
  for (FT x = step; x < limit_x; x += step) {
    for (FT y = y_start; y < limit_y; y += step) {
      const Point_2 query(x, y);

      CGAL::Barycentric_coordinates::mean_value_coordinates_2(
        vertices, query, std::back_inserter(coordinates));
      const FT coordinate_sum =
        coordinates[count + 0] +
        coordinates[count + 1] +
        coordinates[count + 2] +
        coordinates[count + 3] +
        coordinates[count + 4] +
        coordinates[count + 5] +
        coordinates[count + 6] +
        coordinates[count + 7] +
        coordinates[count + 8] +
        coordinates[count + 9] ;
      const Point_2 linear_combination(
        vertices[0].x() * coordinates[count + 0] +
        vertices[1].x() * coordinates[count + 1] +
        vertices[2].x() * coordinates[count + 2] +
        vertices[3].x() * coordinates[count + 3] +
        vertices[4].x() * coordinates[count + 4] +
        vertices[5].x() * coordinates[count + 5] +
        vertices[6].x() * coordinates[count + 6] +
        vertices[7].x() * coordinates[count + 7] +
        vertices[8].x() * coordinates[count + 8] +
        vertices[9].x() * coordinates[count + 9] ,
        vertices[0].y() * coordinates[count + 0] +
        vertices[1].y() * coordinates[count + 1] +
        vertices[2].y() * coordinates[count + 2] +
        vertices[3].y() * coordinates[count + 3] +
        vertices[4].y() * coordinates[count + 4] +
        vertices[5].y() * coordinates[count + 5] +
        vertices[6].y() * coordinates[count + 6] +
        vertices[7].y() * coordinates[count + 7] +
        vertices[8].y() * coordinates[count + 8] +
        vertices[9].y() * coordinates[count + 9] );
      const Point_2 difference(
        linear_combination.x() - query.x(),
        linear_combination.y() - query.y());
      assert(
        CGAL::abs(coordinate_sum - FT(1)) < epsilon &&
        CGAL::abs(difference.x()) < epsilon &&
        CGAL::abs(difference.y()) < epsilon );
      count += 10;
    }
  }
}

int main() {

  test_mv_const_linear_precision< CGAL::Simple_cartesian<double> >();
  test_mv_const_linear_precision< CGAL::Exact_predicates_exact_constructions_kernel >();
  test_mv_const_linear_precision< CGAL::Exact_predicates_inexact_constructions_kernel >();

  std::cout << "test_mv_const_linear_precision: PASSED" << std::endl;
  return EXIT_SUCCESS;
}
