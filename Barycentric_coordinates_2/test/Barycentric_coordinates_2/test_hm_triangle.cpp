#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/Delaunay_domain_2.h>
#include <CGAL/Barycentric_coordinates_2/Harmonic_coordinates_2.h>
#include <CGAL/Barycentric_coordinates_2/triangle_coordinates_2.h>

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
    Point_2(0, 1)
  };

  std::vector<FT> tri_coordinates;
  std::vector<FT>  hm_coordinates;

  std::list<Point_2> seeds;
  seeds.push_back(Point_2(FT(1) / FT(10), FT(1) / FT(10)));

  Domain domain(vertices);
  domain.create(FT(1) / FT(100), seeds);

  HMC2 harmonic_coordinates_2(vertices, domain);
  harmonic_coordinates_2.compute();

  std::size_t count = 0;
  const FT epsilon = FT(1) / FT(100000);
  for (std::size_t i = 0; i < domain.number_of_vertices(); ++i) {
    const auto& query = domain.vertex(i);

    CGAL::Barycentric_coordinates::triangle_coordinates_2(
      vertices[0], vertices[1], vertices[2], query, std::back_inserter(tri_coordinates));
    harmonic_coordinates_2(
      i, std::back_inserter(hm_coordinates));

    assert(hm_coordinates[count + 0] >= FT(0) && hm_coordinates[count + 0] <= FT(1));
    assert(hm_coordinates[count + 1] >= FT(0) && hm_coordinates[count + 1] <= FT(1));
    assert(hm_coordinates[count + 2] >= FT(0) && hm_coordinates[count + 2] <= FT(1));

    assert(
      CGAL::abs(tri_coordinates[count + 0] - hm_coordinates[count + 0]) < epsilon &&
      CGAL::abs(tri_coordinates[count + 1] - hm_coordinates[count + 1]) < epsilon &&
      CGAL::abs(tri_coordinates[count + 2] - hm_coordinates[count + 2]) < epsilon );
    count += 3;
  }

  std::cout << "test_hm_triangle: PASSED" << std::endl;
  return EXIT_SUCCESS;
}
