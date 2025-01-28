#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/Wachspress_coordinates_2.h>
#include <CGAL/Barycentric_coordinates_2/Discrete_harmonic_coordinates_2.h>

using Kernel  = CGAL::Exact_predicates_exact_constructions_kernel;
using FT      = typename Kernel::FT;
using Point_2 = typename Kernel::Point_2;

int main() {

  const std::vector<Point_2> vertices = {
    Point_2(0, 0),
    Point_2(1, 0),
    Point_2(1, 1),
    Point_2(0, 1)
  };

  std::vector<FT> wp_coordinates;
  std::vector<FT> dh_coordinates;

  const FT step  = FT(1) / FT(100);
  const FT scale = FT(100);

  std::size_t count = 0;
  const FT limit = scale * step;
  for (FT x = step; x < limit; x += step) {
    for (FT y = step; y < limit; y += step) {
      const Point_2 query(x, y);

      CGAL::Barycentric_coordinates::wachspress_coordinates_2(
        vertices, query, std::back_inserter(wp_coordinates));
      CGAL::Barycentric_coordinates::discrete_harmonic_coordinates_2(
        vertices, query, std::back_inserter(dh_coordinates));

      assert(wp_coordinates[count + 0] >= FT(0) && wp_coordinates[count + 0] <= FT(1));
      assert(wp_coordinates[count + 1] >= FT(0) && wp_coordinates[count + 1] <= FT(1));
      assert(wp_coordinates[count + 2] >= FT(0) && wp_coordinates[count + 2] <= FT(1));
      assert(wp_coordinates[count + 3] >= FT(0) && wp_coordinates[count + 3] <= FT(1));

      assert(dh_coordinates[count + 0] >= FT(0) && dh_coordinates[count + 0] <= FT(1));
      assert(dh_coordinates[count + 1] >= FT(0) && dh_coordinates[count + 1] <= FT(1));
      assert(dh_coordinates[count + 2] >= FT(0) && dh_coordinates[count + 2] <= FT(1));
      assert(dh_coordinates[count + 3] >= FT(0) && dh_coordinates[count + 3] <= FT(1));

      assert(
        CGAL::abs(wp_coordinates[count + 0] - dh_coordinates[count + 0]) == FT(0) &&
        CGAL::abs(wp_coordinates[count + 1] - dh_coordinates[count + 1]) == FT(0) &&
        CGAL::abs(wp_coordinates[count + 2] - dh_coordinates[count + 2]) == FT(0) &&
        CGAL::abs(wp_coordinates[count + 3] - dh_coordinates[count + 3]) == FT(0) );
      count += 4;
    }
  }

  std::cout << "test_wp_dh_unit_square: PASSED" << std::endl;
  return EXIT_SUCCESS;
}
