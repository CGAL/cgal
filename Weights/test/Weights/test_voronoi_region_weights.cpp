#include <CGAL/Simple_cartesian.h>
#include <CGAL/Weights/voronoi_region_weights.h>

// Typedefs.
using Kernel  = CGAL::Simple_cartesian<double>;
using Point_2 = typename Kernel::Point_2;
using Point_3 = typename Kernel::Point_3;

int main() {

  // 2D configuration.
  const Point_2 p2 = Point_2(-1, 0);
  const Point_2 q2 = Point_2( 0, 0);
  const Point_2 r2 = Point_2( 0, 1);

  // 3D configuration.
  const Point_3 p3 = Point_3(-1, 0, 1);
  const Point_3 q3 = Point_3( 0, 0, 1);
  const Point_3 r3 = Point_3( 0, 1, 1);

  // Compute weights.
  std::cout << "2D voronoi: " <<
    CGAL::Weights::voronoi_area(p2, q2, r2) << std::endl;
  std::cout << "3D voronoi: " <<
    CGAL::Weights::voronoi_area(p3, q3, r3) << std::endl;

  return EXIT_SUCCESS;
}
