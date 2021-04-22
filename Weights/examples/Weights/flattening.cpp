#include <CGAL/Simple_cartesian.h>
#include <CGAL/Weights/authalic_weights.h>
#include <CGAL/Weights/three_point_family_weights.h>

// Typedefs.
using Kernel  = CGAL::Simple_cartesian<double>;
using FT      = typename Kernel::FT;
using Point_3 = typename Kernel::Point_3;

int main() {
  std::cout << std::fixed;

  // 3D configuration.
  const Point_3 p0(0, 1, 1);
  const Point_3 p1(2, 0, 1);
  const Point_3 p2(7, 1, 1);
  const Point_3 q0(3, 1, 1);

  // Choose a type of the weight:
  // e.g. 0 - Wachspress (WP) weight; 1 - mean value (MV);
  const FT wp = FT(0);
  const FT mv = FT(1);

  // Compute WP and MV weights.
  std::cout << "3D Wachspress (WP, q0): ";
  std::cout << CGAL::Weights::three_point_family_weight(p0, p1, p2, q0, wp) << std::endl;
  std::cout << "3D mean value (MV, q0): ";
  std::cout << CGAL::Weights::three_point_family_weight(p0, p1, p2, q0, mv) << std::endl;

  // Converge WP towards MV.
  std::cout << "Converge WP to MV on q0: " << std::endl;
  const FT step = FT(1) / FT(10);
  for (FT x = FT(0); x <= FT(1); x += step) {
    std::cout << "3D x: ";
    std::cout << CGAL::Weights::three_point_family_weight(p0, p1, p2, q0, x) << std::endl;
  }

  // Compute WP weights for q1, which is not on the plane [p0, p1, p2].
  Point_3 q1(3, 1, 2);
  std::cout << "3D wachspress (WP, q1): ";
  std::cout << CGAL::Weights::three_point_family_weight(p0, p1, p2, q1, wp) << std::endl;

  // Converge q1 towards q0 that is we flatten the configuration.
  // We also compare the result with the authalic weight.
  std::cout << "Converge q1 to q0: " << std::endl;
  for (FT x = FT(0); x <= FT(1); x += step) {
    std::cout << "3D wachspress/authalic: ";
    q1 = Point_3(3, 1, FT(2) - x);
    std::cout << CGAL::Weights::three_point_family_weight(p0, p1, p2, q1, wp) << "/";
    std::cout << CGAL::Weights::authalic_weight(p0, p1, p2, q1) << std::endl;
  }
  return EXIT_SUCCESS;
}
