#include <CGAL/Simple_cartesian.h>
#include <CGAL/Weights/inverse_distance_weights.h>

// Typedefs.
using Kernel  = CGAL::Simple_cartesian<double>;
using Point_2 = typename Kernel::Point_2;
using Point_3 = typename Kernel::Point_3;

// Custom traits class that has only two objects.
struct Custom_traits {
  using FT = typename Kernel::FT;
  using Point_2 = typename Kernel::Point_2;
  using Point_3 = typename Kernel::Point_3;

  decltype(auto) compute_squared_distance_2_object() const {
    return Kernel::Compute_squared_distance_2();
  }
  decltype(auto) compute_squared_distance_3_object() const {
    return Kernel::Compute_squared_distance_3();
  }
};

int main() {

  // 2D configuration.
  const Point_2 p2(0, 0);
  const Point_2 q2(0, 1);

  // 3D configuration.
  const Point_3 p3(0, 0, 1);
  const Point_3 q3(0, 1, 1);

  // Create custom traits.
  const Custom_traits ctraits;

  // Compute weights.
  std::cout << "2D/3D weight: ";
  std::cout << CGAL::Weights::inverse_distance_weight(p2, q2, ctraits) << "/";
  std::cout << CGAL::Weights::inverse_distance_weight(p3, q3, ctraits) << std::endl;
  return EXIT_SUCCESS;
}
