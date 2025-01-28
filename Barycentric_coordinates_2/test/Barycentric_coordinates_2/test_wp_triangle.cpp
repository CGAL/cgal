#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/triangle_coordinates_2.h>
#include <CGAL/Barycentric_coordinates_2/Wachspress_coordinates_2.h>

using Kernel  = CGAL::Exact_predicates_exact_constructions_kernel;
using FT      = typename Kernel::FT;
using Point_2 = typename Kernel::Point_2;

int main() {

  const std::vector<Point_2> vertices = {
    Point_2(0, 0),
    Point_2(1, 0),
    Point_2(0, 1)
  };

  std::vector<FT> tri_coordinates;
  std::vector<FT>  wp_coordinates;

  const FT step  = FT(1) / FT(100);
  const FT scale = FT(50);

  std::size_t count = 0;
  const FT limit = scale * step;

  for (FT x = step; x < limit; x += step) {
    for (FT y = step; y < limit; y += step) {
      const Point_2 query(x, y);

      CGAL::Barycentric_coordinates::triangle_coordinates_2(
        vertices[0], vertices[1], vertices[2], query, std::back_inserter(tri_coordinates));
      CGAL::Barycentric_coordinates::wachspress_coordinates_2(
        vertices, query, std::back_inserter(wp_coordinates));

      assert(tri_coordinates[count + 0] >= FT(0) && tri_coordinates[count + 0] <= FT(1));
      assert(tri_coordinates[count + 1] >= FT(0) && tri_coordinates[count + 1] <= FT(1));
      assert(tri_coordinates[count + 2] >= FT(0) && tri_coordinates[count + 2] <= FT(1));

      assert(wp_coordinates[count + 0] >= FT(0) && wp_coordinates[count + 0] <= FT(1));
      assert(wp_coordinates[count + 1] >= FT(0) && wp_coordinates[count + 1] <= FT(1));
      assert(wp_coordinates[count + 2] >= FT(0) && wp_coordinates[count + 2] <= FT(1));

      assert(
        tri_coordinates[count + 0] - wp_coordinates[count + 0] == FT(0) &&
        tri_coordinates[count + 1] - wp_coordinates[count + 1] == FT(0) &&
        tri_coordinates[count + 2] - wp_coordinates[count + 2] == FT(0) );
      count += 3;
    }
  }

  std::cout << "test_wp_triangle: PASSED" << std::endl;
  return EXIT_SUCCESS;
}
