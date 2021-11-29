#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/Mean_value_coordinates_2.h>

// Typedefs.
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;

using FT      = Kernel::FT;
using Point_2 = Kernel::Point_2;
using Policy  = CGAL::Barycentric_coordinates::Computation_policy_2;

int main() {

  // Construct a star-shaped polygon.
  const std::vector<Point_2> star_shaped = {
    Point_2(0.0, 0.0), Point_2( 0.1, -0.8), Point_2(0.3, 0.0), Point_2(0.6, -0.5),
    Point_2(0.6, 0.1), Point_2( 1.1,  0.6), Point_2(0.3, 0.2), Point_2(0.1,  0.8),
    Point_2(0.1, 0.2), Point_2(-0.7,  0.0) };

  // Construct some interior points in the polygon.
  const std::vector<Point_2> interior_points = {
    Point_2(0.12, -0.45), Point_2(0.55, -0.3), Point_2(0.9 , 0.45),
    Point_2(0.15,  0.35), Point_2(-0.4, 0.04), Point_2(0.11, 0.11),
    Point_2(0.28,  0.12), // the only point in the kernel of the star shaped polygon
    Point_2(0.55,  0.11) };

  // Choose a computation policy.
  // We do not check for edge cases since we know
  // that all our points are strictly interior.
  const Policy policy = Policy::PRECISE;

  // Create a vector `std::vector` to store coordinates.
  std::vector<FT> coordinates;
  coordinates.reserve(star_shaped.size());

  // Compute mean value coordinates for all interior points.
  std::cout << std::endl << "mean value coordinates (interior): " << std::endl << std::endl;

  for (const auto& query : interior_points) {
    coordinates.clear();
    CGAL::Barycentric_coordinates::mean_value_coordinates_2(
      star_shaped, query, std::back_inserter(coordinates), policy);

    // Output mean value coordinates.
    for (std::size_t i = 0; i < coordinates.size() - 1; ++i) {
      std::cout << coordinates[i] << ", ";
    }
    std::cout << coordinates[coordinates.size() - 1] << std::endl;
  }
  std::cout << std::endl;

  return EXIT_SUCCESS;
}
