#include <CGAL/Installation/internal/disable_deprecation_warnings_and_errors.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/Generalized_barycentric_coordinates_2.h>
#include <CGAL/Barycentric_coordinates_2/Mean_value_2.h>

// Typedefs.
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;

using FT      = Kernel::FT;
using Point_2 = Kernel::Point_2;

using Mean_value =
  CGAL::Barycentric_coordinates::Mean_value_2<Kernel>;
using Mean_value_coordinates =
  CGAL::Barycentric_coordinates::Generalized_barycentric_coordinates_2<Mean_value, Kernel>;

int main() {

  // Construct a star-shaped polygon.
  const std::vector<Point_2> star_shaped = {
    Point_2(0.0, 0.0), Point_2( 0.1, -0.8), Point_2(0.3, 0.0), Point_2(0.6, -0.5),
    Point_2(0.6, 0.1), Point_2( 1.1,  0.6), Point_2(0.3, 0.2), Point_2(0.1,  0.8),
    Point_2(0.1, 0.2), Point_2(-0.7,  0.0) };

  // Construct the class with mean value coordinates
  // for the star-shaped polygon defined above.
  Mean_value_coordinates mean_value_coordinates(
    star_shaped.begin(), star_shaped.end());

  // Print some information about the polygon and coordinates.
  mean_value_coordinates.print_information();

  // Construct some interior points in the polygon.
  const std::vector<Point_2> interior_points = {
    Point_2(0.12, -0.45), Point_2(0.55, -0.3), Point_2(0.9 , 0.45),
    Point_2(0.15,  0.35), Point_2(-0.4, 0.04), Point_2(0.11, 0.11),
    Point_2(0.28,  0.12), // the only point in the kernel of the star shaped polygon
    Point_2(0.55,  0.11) };

  // We speed up the computation using the O(n) algorithm called with the
  // parameter CGAL::Barycentric_coordinates::FAST.
  // The default one is CGAL::Barycentric_coordinates::PRECISE.
  const auto type_of_algorithm = CGAL::Barycentric_coordinates::FAST;

  // We also speed up the computation by using the parameter
  // query_point_location = CGAL::Barycentric_coordinates::ON_BOUNDED_SIDE.
  const auto query_point_location = CGAL::Barycentric_coordinates::ON_BOUNDED_SIDE;

  // Create a vector `std::vector` to store coordinates.
  std::vector<FT> coordinates;
  coordinates.reserve(star_shaped.size());

  // Compute mean value coordinates for all interior points.
  std::size_t count = 0;
  for (const auto& query : interior_points) {
    coordinates.clear();
    const auto result = mean_value_coordinates(
      query, std::back_inserter(coordinates), query_point_location, type_of_algorithm);

    // Status of the computation.
    const std::string status = (result ? "SUCCESS." : "FAILURE.");
    std::cout << std::endl << "point: " << count << ", status of the computation: " << status << std::endl;
    ++count;

    // Output the coordinates.
    for (std::size_t i = 0; i < coordinates.size() - 1; ++i) {
      std::cout << coordinates[i] << ", ";
    }
    std::cout << coordinates[coordinates.size() - 1] << std::endl;
  }
  std::cout << std::endl;

  return EXIT_SUCCESS;
}
