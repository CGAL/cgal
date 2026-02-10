#include <CGAL/Simple_cartesian.h>
#include <CGAL/Barycentric_coordinates_3/tetrahedron_coordinates.h>

using Kernel = CGAL::Simple_cartesian<double>;

using Point_3 = Kernel::Point_3;

int main() {
  // Construct tetrahedron
  const Point_3 p0(0.0, 0.0, 0.0);
  const Point_3 p1(1.0, 0.0, 0.0);
  const Point_3 p2(0.0, 1.0, 0.0);
  const Point_3 p3(0.0, 0.0, 1.0);

  // Instantiate some interior, boundary, and exterior query points for which we compute coordinates.
  const std::vector<Point_3> queries = {
    Point_3(0.25 , 0.25, 0.25), Point_3(0.3, 0.2, 0.3),         // interior query points
    Point_3(0.1, 0.1, 0.1), Point_3(0.2, 0.5, 0.3),             // interior query points
    Point_3(0.0 , 0.0, 0.5), Point_3(0.4, 0.4, 0.0),            // boundary query points
    Point_3(0.0, 0.4, 0.4), Point_3(0.4, 0.0, 0.4),             // boundary query points
    Point_3(0.5, 0.5, 0.5), Point_3(2.0, 0.0, 0.0),             // exterior query points
    Point_3(-1.0, -1.0, 1.0), Point_3(0.5, 0.5, -2.0)};         // exterior query point

  std::cout << std::endl << "tetrahedron coordinates (all queries): " << std::endl
    << std::endl;

  // Get an array of triangle coordinates for all query points
  for(std::size_t i = 0; i < queries.size(); i++) {
    const auto coords_array =
    CGAL::Barycentric_coordinates::tetrahedron_coordinates(p0, p1, p2, p3, queries[i]);

    std::cout << "tetrahedron coordinates (query " << i << "): " <<
      coords_array[0] << " " << coords_array[1] << " " <<
      coords_array[2] << " " << coords_array[3] << std::endl;
  }

  return EXIT_SUCCESS;
}
