#include <CGAL/Simple_cartesian.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Weights/wachspress_weights.h>

// Typedefs.
using Kernel  = CGAL::Simple_cartesian<double>;
using FT      = typename Kernel::FT;
using Point_2 = typename Kernel::Point_2;

using PointRange = std::vector<Point_2>;
using Creator    = CGAL::Creator_uniform_2<FT, Point_2>;
using Generator  = CGAL::Random_points_in_square_2<Point_2, Creator>;
using Wachspress = CGAL::Weights::Wachspress_weights_2<PointRange, Kernel>;

int main() {

  // Choose how many query points we want to generate.
  const std::size_t num_queries = 10;

  // Create a polygon.
  const std::vector<Point_2> polygon =
    { Point_2(0, 0), Point_2(1, 0), Point_2(1, 1), Point_2(0, 1) };

  // Generate a set of query points.
  std::vector<Point_2> queries;
  queries.reserve(num_queries);
  Generator generator(1.0);
  std::copy_n(generator, num_queries, std::back_inserter(queries));
  assert(queries.size() == num_queries);

  // Allocate memory for weights and coordinates.
  std::vector<FT> weights;
  weights.reserve(num_queries * polygon.size());
  std::vector<FT> coordinates;
  coordinates.reserve(num_queries * polygon.size());

  // Compute barycentric weights.
  Wachspress wachspress(polygon);
  for (const Point_2& query : queries) {
    wachspress(query, std::back_inserter(weights));
  }
  assert(weights.size() == num_queries * polygon.size());

  std::cout << "2D weights: " << std::endl;
  for (std::size_t i = 0; i < weights.size(); i += polygon.size()) {
    for (std::size_t j = 0; j < polygon.size(); ++j) {
      std::cout << weights[i + j] << " ";
    }
    std::cout << std::endl;
  }

  // Normalize weights in order to get barycentric coordinates.
  for (std::size_t i = 0; i < weights.size(); i += polygon.size()) {
    FT sum = 0.0;
    for (std::size_t j = 0; j < polygon.size(); ++j) {
      sum += weights[i + j];
    }
    assert(sum != 0.0);
    for (std::size_t j = 0; j < polygon.size(); ++j) {
      coordinates.push_back(weights[i + j] / sum);
    }
  }
  assert(coordinates.size() == weights.size());

  std::cout << "2D coordinates: " << std::endl;
  for (std::size_t i = 0; i < coordinates.size(); i += polygon.size()) {
    for (std::size_t j = 0; j < polygon.size(); ++j) {
      std::cout << coordinates[i + j] << " ";
    }
    std::cout << std::endl;
  }
  return EXIT_SUCCESS;
}
