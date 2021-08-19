#include <CGAL/Simple_cartesian.h>
#include <CGAL/Barycentric_coordinates_2/triangle_coordinates_2.h>

// Typedefs.
using Kernel = CGAL::Simple_cartesian<double>;

using FT      = Kernel::FT;
using Point_2 = Kernel::Point_2;

int main() {

  // Construct a triangle.
  const Point_2 p0(0.0, 0.0);
  const Point_2 p1(2.0, 0.5);
  const Point_2 p2(1.0, 2.0);

  // Construct several interior, boundary, and exterior query points.
  const std::vector<Point_2> queries = {
    Point_2(0.50, 0.50), // interior query points
    Point_2(1.00, 0.50), Point_2(1.0, 0.75), Point_2(1.00, 1.0),
    Point_2(1.00, 1.25), Point_2(1.0, 1.50), Point_2(0.75, 1.0),
    Point_2(1.25, 1.00), Point_2(1.5, 0.75),

    Point_2(2.0, 0.50), Point_2(1.0, 2.00), // boundary query points
    Point_2(1.0, 0.25), Point_2(1.5, 1.25), Point_2(0.5, 1.0),

    Point_2(0.25, 1.00), Point_2(0.50, 1.75), // exterior query points
    Point_2(1.50, 1.75), Point_2(1.75, 1.50) };

  // Compute triangle coordinates.
  std::vector<FT> coordinates;
  coordinates.reserve(queries.size() * 3);

  for (const auto& query : queries) {
    CGAL::Barycentric_coordinates::triangle_coordinates_2(
      p0, p1, p2, query, std::back_inserter(coordinates));
  }

  // Output all triangle coordinates.
  std::cout << std::endl << "triangle coordinates (all queries): " << std::endl << std::endl;
  for (std::size_t i = 0; i < coordinates.size(); i += 3) {
    std::cout <<
    coordinates[i + 0] << ", " <<
    coordinates[i + 1] << ", " <<
    coordinates[i + 2] << std::endl;
  }
  std::cout << std::endl;
  return EXIT_SUCCESS;
}
