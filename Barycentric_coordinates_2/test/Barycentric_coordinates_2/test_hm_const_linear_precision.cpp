#include <CGAL/Simple_cartesian.h>
#include <CGAL/Barycentric_coordinates_2/Delaunay_domain_2.h>
#include <CGAL/Barycentric_coordinates_2/Harmonic_coordinates_2.h>

using Kernel  = CGAL::Simple_cartesian<double>;
using FT      = typename Kernel::FT;
using Point_2 = typename Kernel::Point_2;

using Vertices = std::vector<Point_2>;
using Domain   = CGAL::Barycentric_coordinates::Delaunay_domain_2<Vertices, Kernel>;
using HMC2     = CGAL::Barycentric_coordinates::Harmonic_coordinates_2<Vertices, Domain, Kernel>;

int main() {

  const Vertices vertices = {
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

  std::vector<Point_2> seeds;
  seeds.push_back(Point_2(0.1, 0.1));

  Domain domain(vertices);
  domain.create(0.01, seeds);

  HMC2 harmonic_coordinates_2(vertices, domain);
  harmonic_coordinates_2.compute();

  std::size_t count = 0;
  std::vector<FT> coordinates;
  const FT epsilon = FT(1) / FT(100000);
  for (std::size_t i = 0; i < domain.number_of_vertices(); ++i) {
    const auto& query = domain.vertex(i);

    harmonic_coordinates_2(
      i, std::back_inserter(coordinates));
    FT coordinate_sum = FT(0);
    for (std::size_t i = 0; i < 50; ++i) {
      assert(
        coordinates[count + i] >= FT(0) &&
        coordinates[count + i] <= FT(1));
      coordinate_sum += coordinates[count + i];
    }
    FT x = FT(0), y = FT(0);
    for (std::size_t i = 0; i < 50; ++i) {
      x += vertices[i].x() * coordinates[count + i];
      y += vertices[i].y() * coordinates[count + i];
    }
    const Point_2 linear_combination(x, y);
    const Point_2 difference(
      linear_combination.x() - query.x(),
      linear_combination.y() - query.y());
    assert(
      CGAL::abs(coordinate_sum - FT(1)) < epsilon &&
      CGAL::abs(difference.x()) < epsilon &&
      CGAL::abs(difference.y()) < epsilon );
    count += 50;
  }

  std::cout << "test_hm_const_linear_precision: PASSED" << std::endl;
  return EXIT_SUCCESS;
}
