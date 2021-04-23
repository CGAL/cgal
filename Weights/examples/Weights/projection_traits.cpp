#include <CGAL/Simple_cartesian.h>
#include <CGAL/Weights/utils.h>
#include <CGAL/Weights/mean_value_weights.h>
#include <CGAL/Weights/internal/Projection_traits_3.h>

// Typedefs.
using Kernel   = CGAL::Simple_cartesian<double>;
using FT       = typename Kernel::FT;
using Point_3  = typename Kernel::Point_3;
using Vector_3 = typename Kernel::Vector_3;
using PTraits  = CGAL::Weights::internal::Projection_traits_3<Kernel>;

int main() {

  // Create a polygon and a query point.
  const std::vector<Point_3> polygon =
    { Point_3(0, 0, 1), Point_3(1, 0, 1), Point_3(1, 1, 0), Point_3(0, 1, 0) };
  const Point_3 query(FT(1) / FT(2), FT(1) / FT(2), FT(1) / FT(2));

  // Allocate memory for weights.
  std::vector<FT> weights;
  weights.reserve(polygon.size());

  // Create projection traits with the right normal.
  const FT quater = FT(1) / FT(4);
  const Vector_3 normal(0, quater, quater);
  const PTraits ptraits(normal);

  // Compute weights.
  CGAL::Weights::mean_value_weights_2(polygon, query, std::back_inserter(weights), ptraits);
  assert(weights.size() == polygon.size());

  std::cout << "2D weights: ";
  for (const FT weight : weights) {
    std::cout << weight << " ";
  }
  std::cout << std::endl;
  return EXIT_SUCCESS;
}
