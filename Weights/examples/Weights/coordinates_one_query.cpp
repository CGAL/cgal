#include <CGAL/Simple_cartesian.h>
#include <CGAL/Weights/discrete_harmonic_weights.h>

// Typedefs.
using Kernel  = CGAL::Simple_cartesian<double>;
using FT      = typename Kernel::FT;
using Point_2 = typename Kernel::Point_2;

int main() {

  // Create a polygon and a query point.
  const std::vector<Point_2> polygon =
    { Point_2(0, 0), Point_2(1, 0), Point_2(1, 1), Point_2(0, 1) };
  const Point_2 query(0.5, 0.5);

  // Allocate memory for weights and coordinates.
  std::vector<FT> weights;
  weights.reserve(polygon.size());
  std::vector<FT> coordinates;
  coordinates.reserve(polygon.size());

  // Compute barycentric weights.
  CGAL::Weights::discrete_harmonic_weights_2(polygon, query, std::back_inserter(weights));
  assert(weights.size() == polygon.size());

  std::cout << "2D weights: ";
  for (const FT weight : weights) {
    std::cout << weight << " ";
  }
  std::cout << std::endl;

  // Normalize weights in order to get barycentric coordinates.
  FT sum = 0.0;
  for (const FT weight : weights) {
    sum += weight;
  }
  assert(sum != 0.0);
  for (const FT weight : weights) {
    coordinates.push_back(weight / sum);
  }
  assert(coordinates.size() == weights.size());

  std::cout << "2D coordinates: ";
  for (const FT coordinate : coordinates) {
    std::cout << coordinate << " ";
  }
  std::cout << std::endl;
  return EXIT_SUCCESS;
}
