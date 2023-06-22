#include <CGAL/Simple_cartesian.h>
#include <CGAL/Weights/three_point_family_weights.h>

// Typedefs.
using Kernel  = CGAL::Simple_cartesian<double>;
using FT      = typename Kernel::FT;
using Point_2 = typename Kernel::Point_2;

int main() {
  std::cout << std::fixed;

  // 3D configuration.
  const Point_2 p0(0, 1);
  const Point_2 p1(2, 0);
  const Point_2 p2(7, 1);
  const Point_2 q0(3, 1);

  // Choose a type of the weight:
  // e.g. 0 - Wachspress (WP) weight; 1 - mean value (MV);
  const FT wp = 0.0;
  const FT mv = 1.0;

  // Compute WP and MV weights.
  std::cout << "3D Wachspress (WP, q0): ";
  std::cout << CGAL::Weights::three_point_family_weight(p0, p1, p2, q0, wp) << std::endl;
  std::cout << "3D mean value (MV, q0): ";
  std::cout << CGAL::Weights::three_point_family_weight(p0, p1, p2, q0, mv) << std::endl;

  // Converge WP towards MV.
  std::cout << "Converge WP to MV on q0: " << std::endl;
  const FT step = 0.1;
  for (FT x = 0.0; x <= 1.0; x += step) {
    std::cout << "3D x: ";
    std::cout << CGAL::Weights::three_point_family_weight(p0, p1, p2, q0, x) << std::endl;
  }
  return EXIT_SUCCESS;
}
