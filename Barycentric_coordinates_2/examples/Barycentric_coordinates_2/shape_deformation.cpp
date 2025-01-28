#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/Delaunay_domain_2.h>
#include <CGAL/Barycentric_coordinates_2/Harmonic_coordinates_2.h>

// Typedefs.
using Kernel      = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT          = Kernel::FT;
using Point_2     = Kernel::Point_2;
using Point_range = std::vector<Point_2>;

using Domain =
  CGAL::Barycentric_coordinates::Delaunay_domain_2<Point_range, Kernel>;
using Harmonic_coordinates_2 =
  CGAL::Barycentric_coordinates::Harmonic_coordinates_2<Point_range, Domain, Kernel>;

int main() {

  // Construct the source and target shapes.
  // The number of vertices in both shapes must be equal.
  const std::vector<Point_2> source_shape = {
    Point_2(1, 0), Point_2(2, 0), Point_2(3, 3), Point_2(4, 0), Point_2(5, 0),
    Point_2(4, 3), Point_2(4, 5), Point_2(5, 4), Point_2(5, 5), Point_2(4, 6),
    Point_2(2, 6), Point_2(1, 5), Point_2(1, 4), Point_2(2, 5), Point_2(2, 3)
  };
  const std::vector<Point_2> target_shape = {
    Point_2(2, 0), Point_2(3, 0), Point_2(3, 3), Point_2(3, 0), Point_2(4, 0),
    Point_2(4, 3), Point_2(4, 5), Point_2(5, 6), Point_2(5, 7), Point_2(4, 6),
    Point_2(2, 6), Point_2(1, 7), Point_2(1, 6), Point_2(2, 5), Point_2(2, 3)
  };
  assert(target_shape.size() == source_shape.size());

  // Use seeds to mark the interior part of the source shape.
  const std::vector<Point_2> seeds = { Point_2(3, 5) };

  // Construct a Delaunay domain.
  const FT max_edge_length = FT(1) / FT(3);
  Domain domain(source_shape);
  domain.create(max_edge_length, seeds);

  // Use it to store coordinates.
  std::vector< std::vector<FT> > coordinates;
  coordinates.reserve(domain.number_of_vertices());

  // Compute harmonic coordinates at the vertices of the
  // discretized interior domain of the source shape.
  Harmonic_coordinates_2 harmonic_coordinates_2(source_shape, domain);
  harmonic_coordinates_2.compute();
  harmonic_coordinates_2(std::back_inserter(coordinates));

  // Deform the source domain into the target domain.
  // We output only the first 20 results.
  std::cout << std::endl << "shape deformation: " << std::endl << std::endl;

  for (std::size_t k = 0; k < 20; ++k) {
    FT x = FT(0), y = FT(0);
    for (std::size_t i = 0; i < coordinates[k].size(); ++i) {
      x += coordinates[k][i] * target_shape[i].x();
      y += coordinates[k][i] * target_shape[i].y();
    }
    std::cout << "deformed domain vertex: (" << x << ", " << y << ")" << std::endl;
  }
  std::cout << std::endl;

  return EXIT_SUCCESS;
}
