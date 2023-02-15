#include <CGAL/Simple_cartesian.h>
#include <CGAL/Barycentric_coordinates_2/Delaunay_domain_2.h>
#include <CGAL/Barycentric_coordinates_2/Harmonic_coordinates_2.h>

// Typedefs.
using Kernel      = CGAL::Simple_cartesian<double>;
using Point_2     = Kernel::Point_2;
using Point_range = std::vector<Point_2>;

using Domain =
  CGAL::Barycentric_coordinates::Delaunay_domain_2<Point_range, Kernel>;
using Harmonic_coordinates_2 =
  CGAL::Barycentric_coordinates::Harmonic_coordinates_2<Point_range, Domain, Kernel>;

int main() {

  // Construct a simple polygon.
  const std::vector<Point_2> polygon = {
    Point_2(0.03, 0.05), Point_2(0.07, 0.04), Point_2(0.10, 0.04),
    Point_2(0.14, 0.04), Point_2(0.17, 0.07), Point_2(0.20, 0.09),
    Point_2(0.22, 0.11), Point_2(0.25, 0.11), Point_2(0.27, 0.10),
    Point_2(0.30, 0.07), Point_2(0.31, 0.04), Point_2(0.34, 0.03),
    Point_2(0.37, 0.02), Point_2(0.40, 0.03), Point_2(0.42, 0.04),
    Point_2(0.44, 0.07), Point_2(0.45, 0.10), Point_2(0.46, 0.13),
    Point_2(0.46, 0.19), Point_2(0.47, 0.26), Point_2(0.47, 0.31),
    Point_2(0.47, 0.35), Point_2(0.45, 0.37), Point_2(0.41, 0.38),
    Point_2(0.38, 0.37), Point_2(0.35, 0.36), Point_2(0.32, 0.35),
    Point_2(0.30, 0.37), Point_2(0.28, 0.39), Point_2(0.25, 0.40),
    Point_2(0.23, 0.39), Point_2(0.21, 0.37), Point_2(0.21, 0.34),
    Point_2(0.23, 0.32), Point_2(0.24, 0.29), Point_2(0.27, 0.24),
    Point_2(0.29, 0.21), Point_2(0.29, 0.18), Point_2(0.26, 0.16),
    Point_2(0.24, 0.17), Point_2(0.23, 0.19), Point_2(0.24, 0.22),
    Point_2(0.24, 0.25), Point_2(0.21, 0.26), Point_2(0.17, 0.26),
    Point_2(0.12, 0.24), Point_2(0.07, 0.20), Point_2(0.03, 0.15),
    Point_2(0.01, 0.10), Point_2(0.02, 0.07)
  };

  // Use seeds to mark the interior part of the polygon.
  std::list<Point_2> seeds;
  seeds.push_back(Point_2(0.1, 0.1));

  // Construct a Delaunay domain.
  const double max_edge_length = 0.01;
  Domain domain(polygon);
  domain.create(max_edge_length, seeds);

  // Compute harmonic coordinates at the vertices of the domain.
  Harmonic_coordinates_2 harmonic_coordinates_2(polygon, domain);
  harmonic_coordinates_2.compute();

  // Use it to store coordinates.
  std::vector<double> coordinates;
  coordinates.reserve(polygon.size());

  // Output harmonic coordinates.
  // We output only the first 20 results.
  std::cout.precision(1);
  std::cout << std::endl << "harmonic coordinates (computed): " << std::endl << std::endl;

  for (std::size_t k = 0; k < 20; ++k) {
    coordinates.clear();
    harmonic_coordinates_2(k, std::back_inserter(coordinates));
    for (std::size_t i = 0; i < coordinates.size() - 1; ++i) {
      std::cout << coordinates[i] << ", ";
    }
    std::cout << coordinates[coordinates.size() - 1] << std::endl;
  }

  // Evaluate harmonic coordinates at the barycenters of the domain triangles.
  // We output only the first 20 results.
  std::cout << std::endl << "harmonic coordinates (evaluated): " << std::endl << std::endl;

  std::vector<Point_2> barycenters;
  domain.barycenters(std::back_inserter(barycenters));

  for (std::size_t k = 0; k < 20; ++k) {
    coordinates.clear();
    harmonic_coordinates_2(barycenters[k], std::back_inserter(coordinates));
    for (std::size_t i = 0; i < coordinates.size() - 1; ++i) {
      std::cout << coordinates[i] << ", ";
    }
    std::cout << coordinates[coordinates.size() - 1] << std::endl;
  }
  std::cout << std::endl;

  return EXIT_SUCCESS;
}
