#include <CGAL/Simple_cartesian.h>
#include <CGAL/Weights/utils.h>
#include <CGAL/Weights/authalic_weights.h>

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
  std::cout << "2D authalic: " <<
    CGAL::Weights::authalic_weight(t2, r2, p2, q2) << std::endl;
  std::cout << "3D authalic: " <<
    CGAL::Weights::authalic_weight(t3, r3, p3, q3) << std::endl;
  std::cout << "--------------" << std::endl;

  // Construct a 2D weight.
  const auto w2 =
    CGAL::Weights::half_authalic_weight(
      CGAL::Weights::cotangent(t2, r2, q2),
      CGAL::Weights::squared_distance(q2, r2)) +
    CGAL::Weights::half_authalic_weight(
      CGAL::Weights::cotangent(q2, r2, p2),
      CGAL::Weights::squared_distance(q2, r2));
  std::cout << "2D authalic: " << w2 << std::endl;

  // Construct a 3D weight.
  const auto w3 =
    CGAL::Weights::half_authalic_weight(
      CGAL::Weights::cotangent(t3, r3, q3),
      CGAL::Weights::squared_distance(q3, r3)) +
    CGAL::Weights::half_authalic_weight(
      CGAL::Weights::cotangent(q3, r3, p3),
      CGAL::Weights::squared_distance(q3, r3));
  std::cout << "3D authalic: " << w3 << std::endl;

  return EXIT_SUCCESS;
}
