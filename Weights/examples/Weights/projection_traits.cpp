#include <CGAL/Simple_cartesian.h>
#include <CGAL/Weights/wachspress_weights.h>
#include <CGAL/Weights/mean_value_weights.h>
#include <CGAL/Projection_traits_xy_3.h>

// Typedefs.
using Kernel   = CGAL::Simple_cartesian<double>;
using FT       = typename Kernel::FT;
using Point_3  = typename Kernel::Point_3;
using Vector_3 = typename Kernel::Vector_3;
using PTraits  = CGAL::Projection_traits_xy_3<Kernel>;

int main() {

  // Create a polygon and a query point.
  const std::vector<Point_3> polygon =
    { Point_3(0, 0, 1), Point_3(1, 0, 1), Point_3(1, 1, 1), Point_3(0, 1, 1) };
  const Point_3 query(0.5, 0.5, 1);

  // Allocate memory for weights.
  std::vector<FT> weights;
  weights.reserve(polygon.size());

  // Create projection traits.
  const PTraits ptraits;

  // Compute weights.
  CGAL::Weights::mean_value_weights_2(polygon, query, std::back_inserter(weights), ptraits);
  assert(weights.size() == polygon.size());

  std::cout << "2D mean value weights: ";
  for (const FT weight : weights) {
    std::cout << weight << " ";
  }
  std::cout << std::endl;

  // Almost coplanar case.
  const FT eps = 0.001;
  const Point_3 t3(-1,  0, 1.0 + eps);
  const Point_3 r3( 0, -1, 1);
  const Point_3 p3( 1,  0, 1.0 + eps);
  const Point_3 q3( 0,  0, 1);
  std::cout << "2D Wachspress weight: " <<
    CGAL::Weights::wachspress_weight(t3, r3, p3, q3, ptraits) << std::endl;
  return EXIT_SUCCESS;
}
