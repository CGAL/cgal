#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/Wachspress_coordinates_2.h>

using Kernel  = CGAL::Exact_predicates_exact_constructions_kernel;
using FT      = typename Kernel::FT;
using Point_2 = typename Kernel::Point_2;

int main() {

  const std::vector<Point_2> vertices = {
    Point_2(0, 0),
    Point_2(1, 0),
    Point_2( FT(7) / FT(4), FT(3) / FT(4)),
    Point_2( FT(5) / FT(4), FT(3) / FT(2)),
    Point_2( FT(1) / FT(4), FT(3) / FT(2)),
    Point_2(-FT(1) / FT(2), FT(5) / FT(4))
  };

  const FT step    = FT(1) / FT(100);
  const FT x_scale = FT(100);
  const FT y_scale = FT(115);

  std::size_t count = 0;
  const Point_2 zero(0, 0);
  const FT limit_x = x_scale * step;
  const FT limit_y = y_scale * step;

  std::vector<FT> coordinates;
  for (FT x = step; x < limit_x; x += step) {
    for (FT y = step; y < limit_y; y += step) {
      const Point_2 query(x, y);

      CGAL::Barycentric_coordinates::wachspress_coordinates_2(
        vertices, query, std::back_inserter(coordinates));
      const FT coordinate_sum =
        coordinates[count + 0] +
        coordinates[count + 1] +
        coordinates[count + 2] +
        coordinates[count + 3] +
        coordinates[count + 4] +
        coordinates[count + 5] ;
      const Point_2 linear_combination(
        vertices[0].x() * coordinates[count + 0] +
        vertices[1].x() * coordinates[count + 1] +
        vertices[2].x() * coordinates[count + 2] +
        vertices[3].x() * coordinates[count + 3] +
        vertices[4].x() * coordinates[count + 4] +
        vertices[5].x() * coordinates[count + 5] ,
        vertices[0].y() * coordinates[count + 0] +
        vertices[1].y() * coordinates[count + 1] +
        vertices[2].y() * coordinates[count + 2] +
        vertices[3].y() * coordinates[count + 3] +
        vertices[4].y() * coordinates[count + 4] +
        vertices[5].y() * coordinates[count + 5] );
      const Point_2 difference(
        linear_combination.x() - query.x(),
        linear_combination.y() - query.y());
      assert( (coordinate_sum == FT(1)) && (difference == zero) );
      count += 6;
    }
  }

  std::cout << "test_wp_const_linear_precision: PASSED" << std::endl;
  return EXIT_SUCCESS;
}
