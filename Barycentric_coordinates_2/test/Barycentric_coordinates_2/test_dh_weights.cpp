#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/Discrete_harmonic_coordinates_2.h>

using Kernel  = CGAL::Exact_predicates_exact_constructions_kernel;
using FT      = typename Kernel::FT;
using Point_2 = typename Kernel::Point_2;

int main() {

  const std::vector<Point_2> vertices = {
    Point_2(0, 0),
    Point_2(1, 0),
    Point_2( FT(5) / FT(4), FT(3) / FT(4)),
    Point_2( FT(1) / FT(2), FT(3) / FT(2)),
    Point_2(-FT(1) / FT(4), FT(3) / FT(4))
  };

  std::vector<FT> weights;
  std::vector<FT> coordinates;
  std::vector<FT> expected;

  const FT step  = FT(1) / FT(100);
  const FT scale = FT(100);

  std::size_t count = 0;
  const FT limit = scale * step;

  for (FT x = step; x < limit; x += step) {
    for (FT y = step; y < limit; y += step) {
      const Point_2 query(x, y);

      CGAL::Barycentric_coordinates::discrete_harmonic_weights_2(
        vertices, query, std::back_inserter(weights));

      FT W = FT(0);
      for (std::size_t j = 0; j < 5; ++j)
        W += weights[count + j];
      const FT inv_W = FT(1) / W;
      for (std::size_t j = 0; j < 5; ++j)
        coordinates.push_back(weights[count + j] * inv_W);

      CGAL::Barycentric_coordinates::discrete_harmonic_coordinates_2(
        vertices, query, std::back_inserter(expected));
      assert(
        coordinates[count + 0] - expected[count + 0] == FT(0) &&
        coordinates[count + 1] - expected[count + 1] == FT(0) &&
        coordinates[count + 2] - expected[count + 2] == FT(0) &&
        coordinates[count + 3] - expected[count + 3] == FT(0) &&
        coordinates[count + 4] - expected[count + 4] == FT(0) );
      count += 5;
    }
  }

  const Point_2 queries[5] = {
    Point_2(FT(1) / FT(std::pow(10.0, 300.0)), FT(1) / FT(std::pow(10.0, 300.0))),
    Point_2(FT(1) / FT(4), FT(1) / FT(std::pow(10.0, 300.0))),
    Point_2(FT(1) / FT(2), FT(1) / FT(std::pow(10.0, 300.0))),
    Point_2(FT(3) / FT(4), FT(1) / FT(std::pow(10.0, 300.0))),
    Point_2(FT(1) - FT(1) / FT(std::pow(10.0, 300.0)), FT(1) / FT(std::pow(10.0, 300.0)))
  };

  for (std::size_t i = 0; i < 5; ++i) {
    CGAL::Barycentric_coordinates::discrete_harmonic_weights_2(
      vertices, queries[i], std::back_inserter(weights));

    FT W = FT(0);
    for(std::size_t j = 0; j < 5; ++j)
      W += weights[count + j];
    const FT inv_W = FT(1) / W;
    for (std::size_t j = 0; j < 5; ++j)
      coordinates.push_back(weights[count + j] * inv_W);

    CGAL::Barycentric_coordinates::discrete_harmonic_coordinates_2(
      vertices, queries[i], std::back_inserter(expected));
    assert(
      coordinates[count + 0] - expected[count + 0] == FT(0) &&
      coordinates[count + 1] - expected[count + 1] == FT(0) &&
      coordinates[count + 2] - expected[count + 2] == FT(0) &&
      coordinates[count + 3] - expected[count + 3] == FT(0) &&
      coordinates[count + 4] - expected[count + 4] == FT(0) );
    count += 5;
  }

  std::cout << "test_dh_weights: PASSED" << std::endl;
  return EXIT_SUCCESS;
}
