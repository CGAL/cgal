#include <CGAL/Simple_cartesian.h>
#include <CGAL/Weights/utils.h>
#include <CGAL/Weights/tangent_weights.h>

// Typedefs.
using Kernel  = CGAL::Simple_cartesian<double>;
using Point_2 = typename Kernel::Point_2;
using Point_3 = typename Kernel::Point_3;

int main() {

  // 2D configuration.
  const Point_2 t2 = Point_2(-1,  0);
  const Point_2 r2 = Point_2( 0, -1);
  const Point_2 p2 = Point_2( 1,  0);
  const Point_2 q2 = Point_2( 0,  0);

  // 3D configuration.
  const Point_3 t3 = Point_3(-1,  0, 1);
  const Point_3 r3 = Point_3( 0, -1, 1);
  const Point_3 p3 = Point_3( 1,  0, 1);
  const Point_3 q3 = Point_3( 0,  0, 1);

  // Compute weights.
  std::cout << "2D tangent: " <<
    CGAL::Weights::tangent_weight(t2, r2, p2, q2) << std::endl;
  std::cout << "3D tangent: " <<
    CGAL::Weights::tangent_weight(t3, r3, p3, q3) << std::endl;
  std::cout << "-------------" << std::endl;

  // Construct a 2D weight.
  const auto w2 =
    CGAL::Weights::half_tangent_weight(
      CGAL::Weights::distance(r2, q2),
      CGAL::Weights::distance(t2, q2),
      CGAL::Weights::area(r2, q2, t2),
      CGAL::Weights::scalar_product(r2, q2, t2)) +
    CGAL::Weights::half_tangent_weight(
      CGAL::Weights::distance(r2, q2),
      CGAL::Weights::distance(p2, q2),
      CGAL::Weights::area(p2, q2, r2),
      CGAL::Weights::scalar_product(p2, q2, r2));
  std::cout << "2D tangent: " << w2 << std::endl;

  // Construct a 3D weight.
  const auto w3 =
    CGAL::Weights::half_tangent_weight(
      CGAL::Weights::tangent_half_angle(
        CGAL::Weights::distance(r3, q3),
        CGAL::Weights::distance(t3, q3),
        CGAL::Weights::area(r3, q3, t3),
        CGAL::Weights::scalar_product(r3, q3, t3)),
      CGAL::Weights::distance(r3, q3)) +
    CGAL::Weights::half_tangent_weight(
      CGAL::Weights::tangent_half_angle(
        CGAL::Weights::distance(r3, q3),
        CGAL::Weights::distance(p3, q3),
        CGAL::Weights::area(p3, q3, r3),
        CGAL::Weights::scalar_product(p3, q3, r3)),
      CGAL::Weights::distance(r3, q3));
  std::cout << "3D tangent: " << w3 << std::endl;

  return EXIT_SUCCESS;
}
