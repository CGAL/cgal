#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/interpolation_functions.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/Delaunay_domain_2.h>
#include <CGAL/Barycentric_coordinates_2/Mean_value_coordinates_2.h>

// Typedefs.
using Kernel     = CGAL::Exact_predicates_inexact_constructions_kernel;
using Projection = CGAL::Projection_traits_xy_3<Kernel>;

using FT          = typename Projection::FT;
using Point       = typename Projection::Point_2;
using Point_range = std::vector<Point>;

using Domain =
  CGAL::Barycentric_coordinates::Delaunay_domain_2<Point_range, Projection>;
using Mean_value_coordinates_2 =
  CGAL::Barycentric_coordinates::Mean_value_coordinates_2<Point_range, Projection>;

using Vertex_function_value = std::map<Point, FT, typename Projection::Less_xy_2>;
using Function_value_access = CGAL::Data_access<Vertex_function_value>;
using Point_with_coordinate = std::pair<Point, FT>;

int main() {

  // Construct a polygon that bounds a three-dimensional terrain.
  // Note that the z-coordinate of each vertex represents the height function.
  // Projection in 2D is performed automatically by the Projection traits class.
  const std::vector<Point> polygon = {
    Point(0.03, 0.05, 0.00), Point(0.07, 0.04, 0.02), Point(0.10, 0.04, 0.04),
    Point(0.14, 0.04, 0.06), Point(0.17, 0.07, 0.08), Point(0.20, 0.09, 0.10),
    Point(0.22, 0.11, 0.12), Point(0.25, 0.11, 0.14), Point(0.27, 0.10, 0.16),
    Point(0.30, 0.07, 0.18), Point(0.31, 0.04, 0.20), Point(0.34, 0.03, 0.22),
    Point(0.37, 0.02, 0.24), Point(0.40, 0.03, 0.26), Point(0.42, 0.04, 0.28),
    Point(0.44, 0.07, 0.30), Point(0.45, 0.10, 0.32), Point(0.46, 0.13, 0.34),
    Point(0.46, 0.19, 0.36), Point(0.47, 0.26, 0.38), Point(0.47, 0.31, 0.40),
    Point(0.47, 0.35, 0.42), Point(0.45, 0.37, 0.44), Point(0.41, 0.38, 0.46),
    Point(0.38, 0.37, 0.48), Point(0.35, 0.36, 0.50), Point(0.32, 0.35, 0.52),
    Point(0.30, 0.37, 0.54), Point(0.28, 0.39, 0.56), Point(0.25, 0.40, 0.58),
    Point(0.23, 0.39, 0.60), Point(0.21, 0.37, 0.62), Point(0.21, 0.34, 0.64),
    Point(0.23, 0.32, 0.66), Point(0.24, 0.29, 0.68), Point(0.27, 0.24, 0.70),
    Point(0.29, 0.21, 0.72), Point(0.29, 0.18, 0.74), Point(0.26, 0.16, 0.76),
    Point(0.24, 0.17, 0.78), Point(0.23, 0.19, 0.80), Point(0.24, 0.22, 0.82),
    Point(0.24, 0.25, 0.84), Point(0.21, 0.26, 0.86), Point(0.17, 0.26, 0.88),
    Point(0.12, 0.24, 0.90), Point(0.07, 0.20, 0.92), Point(0.03, 0.15, 0.94),
    Point(0.01, 0.10, 0.97), Point(0.02, 0.07, 1.00)
  };

  // Construct a Delaunay domain.
  std::vector<Point> seeds;
  seeds.push_back(Point(0.1, 0.1, 0.0));

  Domain domain(polygon);
  domain.create(0.05, seeds);

  // Associate each polygon vertex with the corresponding function value.
  Vertex_function_value vertex_function_value;
  for (const auto& vertex : polygon) {
    vertex_function_value.insert(
      std::make_pair(vertex, vertex.z()));
  }

  // Construct the class with the mean value weights.
  Mean_value_coordinates_2 mean_value_coordinates_2(polygon);

  // Compute mean value coordinates and use them to interpolate data
  // from the polygon boundary to its interior.
  std::vector<FT> coordinates;
  coordinates.reserve(polygon.size());

  std::vector<Point_with_coordinate> boundary;
  boundary.resize(polygon.size());

  std::vector<Point> queries;
  queries.reserve(domain.number_of_vertices());

  for (std::size_t i = 0; i < domain.number_of_vertices(); ++i) {
    const auto& query = domain.vertex(i);

    coordinates.clear();
    mean_value_coordinates_2(query, std::back_inserter(coordinates));
    for (std::size_t i = 0; i < polygon.size(); ++i) {
      boundary[i] = std::make_pair(polygon[i], coordinates[i]);
    }

    const FT f = CGAL::linear_interpolation(
      boundary.begin(), boundary.end(), FT(1),
      Function_value_access(vertex_function_value));
    queries.push_back(Point(query.x(), query.y(), f));
  }

  // Output interpolated heights.
  std::cout << std::endl << "interpolated heights (all queries): " << std::endl << std::endl;
  for (const auto& query : queries) {
    std::cout << query.z() << std::endl;
  }
  std::cout << std::endl;

  return EXIT_SUCCESS;
}
