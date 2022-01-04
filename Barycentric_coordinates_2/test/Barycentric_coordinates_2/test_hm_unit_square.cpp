#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/Delaunay_domain_2.h>
#include <CGAL/Barycentric_coordinates_2/Harmonic_coordinates_2.h>

using Kernel  = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT      = typename Kernel::FT;
using Point_2 = typename Kernel::Point_2;

using Vertices = std::vector<Point_2>;
using Domain   = CGAL::Barycentric_coordinates::Delaunay_domain_2<Vertices, Kernel>;
using HMC2     = CGAL::Barycentric_coordinates::Harmonic_coordinates_2<Vertices, Domain, Kernel>;

int main() {

  const Vertices vertices = {
    Point_2(0, 0),
    Point_2(1, 0),
    Point_2(1, 1),
    Point_2(0, 1)
  };

  std::list<Point_2> seeds;
  seeds.push_back(Point_2(FT(1) / FT(2), FT(1) / FT(2)));

  Domain domain(vertices);
  domain.release_memory();
  domain.create(FT(1) / FT(100), seeds);

  HMC2 harmonic_coordinates_2(vertices, domain);
  harmonic_coordinates_2.clear();
  harmonic_coordinates_2.setup();
  harmonic_coordinates_2.factorize();
  harmonic_coordinates_2.solve();

  std::vector<FT> coordinates;
  const FT epsilon = FT(1) / FT(100000);
  for (std::size_t i = 0; i < domain.number_of_vertices(); ++i) {
    const auto& query = domain.vertex(i);
    const FT ref_value = query.x() * query.y();

    coordinates.clear();
    harmonic_coordinates_2(i, std::back_inserter(coordinates));
    assert(coordinates[0] >= FT(0) && coordinates[0] <= FT(1));
    assert(coordinates[1] >= FT(0) && coordinates[1] <= FT(1));
    assert(coordinates[2] >= FT(0) && coordinates[2] <= FT(1));
    assert(coordinates[3] >= FT(0) && coordinates[3] <= FT(1));
    const FT value = coordinates[2];
    assert(CGAL::abs(ref_value - value) < epsilon);
  }

  harmonic_coordinates_2.release_memory();
  harmonic_coordinates_2.compute();

  std::vector<Point_2> barycenters;
  domain.barycenters(std::back_inserter(barycenters));
  for (const auto& barycenter : barycenters) {
    const FT ref_value = barycenter.x() * barycenter.y();

    coordinates.clear();
    harmonic_coordinates_2(barycenter, std::back_inserter(coordinates));
    assert(coordinates[0] >= FT(0) && coordinates[0] <= FT(1));
    assert(coordinates[1] >= FT(0) && coordinates[1] <= FT(1));
    assert(coordinates[2] >= FT(0) && coordinates[2] <= FT(1));
    assert(coordinates[3] >= FT(0) && coordinates[3] <= FT(1));
    const FT value = coordinates[2];
    assert(CGAL::abs(ref_value - value) < epsilon);
  }

  const FT quater = FT(1) / FT(4);
  coordinates.clear();
  const Point_2 query(FT(1) / FT(2), FT(1) / FT(2));
  harmonic_coordinates_2(query, std::back_inserter(coordinates));
  assert(
    CGAL::abs(coordinates[0] - quater) < epsilon &&
    CGAL::abs(coordinates[1] - quater) < epsilon &&
    CGAL::abs(coordinates[2] - quater) < epsilon &&
    CGAL::abs(coordinates[3] - quater) < epsilon );

  std::cout << "test_hm_unit_square: PASSED" << std::endl;
  return EXIT_SUCCESS;
}
