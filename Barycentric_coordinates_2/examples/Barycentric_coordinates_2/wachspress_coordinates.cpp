#include <CGAL/convex_hull_2.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Barycentric_coordinates_2/Wachspress_coordinates_2.h>

// Typedefs.
using Kernel = CGAL::Simple_cartesian<double>;

using FT      = Kernel::FT;
using Point_2 = Kernel::Point_2;

using Creator   = CGAL::Creator_uniform_2<FT, Point_2>;
using Generator = CGAL::Random_points_in_square_2<Point_2, Creator>;

int main() {

  // Choose how many query points we want to generate.
  const std::size_t num_queries = 100;

  // Create vectors to store query points and polygon vertices.
  std::vector<Point_2> queries, convex;

  // Generate a set of query points.
  queries.reserve(num_queries);
  Generator generator(1.0);
  std::copy_n(generator, num_queries, std::back_inserter(queries));

  // Find the convex hull of the generated query points.
  // This convex hull gives the vertices of a convex polygon
  // that contains all the generated points.
  CGAL::convex_hull_2(
    queries.begin(), queries.end(), std::back_inserter(convex));

  // Compute Wachspress coordinates for all query points.
  std::cout << std::endl << "Wachspress coordinates (interior + boundary): " << std::endl << std::endl;

  std::vector<FT> coordinates;
  coordinates.reserve(convex.size());

  for (const auto& query : queries) {
    coordinates.clear();
    CGAL::Barycentric_coordinates::wachspress_coordinates_2(
      convex, query, std::back_inserter(coordinates));

    for (std::size_t i = 0; i < coordinates.size() - 1; ++i) {
      std::cout << coordinates[i] << ", ";
    }
    std::cout << coordinates[coordinates.size() - 1] << std::endl;
  }
  std::cout << std::endl;

  return EXIT_SUCCESS;
}
