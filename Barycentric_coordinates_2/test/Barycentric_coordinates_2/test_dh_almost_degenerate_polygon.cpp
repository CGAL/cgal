#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/Discrete_harmonic_coordinates_2.h>

using Kernel  = CGAL::Exact_predicates_exact_constructions_kernel;
using FT      = typename Kernel::FT;
using Point_2 = typename Kernel::Point_2;

int main() {

  const std::vector<Point_2> vertices = {
    Point_2(0, 0),
    Point_2(1, 0),
    Point_2(1, FT(1) / FT(std::pow(10.0, 200.0))),
    Point_2(0, FT(1) / FT(std::pow(10.0, 200.0)))
  };

  const Point_2 queries[31] = {
    Point_2(FT(1) / FT(std::pow(10.0, 300.0)), FT(1) / FT(std::pow(10.0, 300.0))) ,
    Point_2(FT(1) / FT(std::pow(10.0, 300.0)), FT(1) / FT(std::pow(10.0, 280.0))) ,
    Point_2(FT(1) / FT(std::pow(10.0, 300.0)), FT(1) / FT(std::pow(10.0, 260.0))) ,
    Point_2(FT(1) / FT(std::pow(10.0, 300.0)), FT(1) / FT(std::pow(10.0, 240.0))) ,
    Point_2(FT(1) / FT(std::pow(10.0, 300.0)), FT(1) / FT(std::pow(10.0, 220.0))) ,
    Point_2(FT(1) / FT(std::pow(10.0, 300.0)), FT(1) / FT(std::pow(10.0, 210.0))) ,
    Point_2(FT(1) / FT(std::pow(10.0, 300.0)), FT(1) / FT(std::pow(10.0, 200.8))) ,
    Point_2(FT(1) / FT(std::pow(10.0, 300.0)), FT(1) / FT(std::pow(10.0, 200.5))) ,
    Point_2(FT(1) / FT(std::pow(10.0, 300.0)), FT(1) / FT(std::pow(10.0, 200.2))) ,
    Point_2(FT(1) / FT(std::pow(10.0, 300.0)), FT(1) / FT(std::pow(10.0, 200.05))),

    Point_2(FT(1) / FT(2), FT(1) / FT(std::pow(10.0, 300.0))) ,
    Point_2(FT(1) / FT(2), FT(1) / FT(std::pow(10.0, 280.0))) ,
    Point_2(FT(1) / FT(2), FT(1) / FT(std::pow(10.0, 260.0))) ,
    Point_2(FT(1) / FT(2), FT(1) / FT(std::pow(10.0, 240.0))) ,
    Point_2(FT(1) / FT(2), FT(1) / FT(std::pow(10.0, 220.0))) ,
    Point_2(FT(1) / FT(2), FT(1) / FT(std::pow(10.0, 210.0))) ,
    Point_2(FT(1) / FT(2), FT(1) / FT(std::pow(10.0, 200.8))) ,
    Point_2(FT(1) / FT(2), FT(1) / FT(std::pow(10.0, 200.5))) ,
    Point_2(FT(1) / FT(2), FT(1) / FT(std::pow(10.0, 200.2))) ,
    Point_2(FT(1) / FT(2), FT(1) / FT(std::pow(10.0, 200.05))),

    Point_2(FT(1) - FT(1) / FT(std::pow(10.0, 300.0)), FT(1) / FT(std::pow(10.0, 300.0))) ,
    Point_2(FT(1) - FT(1) / FT(std::pow(10.0, 300.0)), FT(1) / FT(std::pow(10.0, 280.0))) ,
    Point_2(FT(1) - FT(1) / FT(std::pow(10.0, 300.0)), FT(1) / FT(std::pow(10.0, 260.0))) ,
    Point_2(FT(1) - FT(1) / FT(std::pow(10.0, 300.0)), FT(1) / FT(std::pow(10.0, 240.0))) ,
    Point_2(FT(1) - FT(1) / FT(std::pow(10.0, 300.0)), FT(1) / FT(std::pow(10.0, 220.0))) ,
    Point_2(FT(1) - FT(1) / FT(std::pow(10.0, 300.0)), FT(1) / FT(std::pow(10.0, 210.0))) ,
    Point_2(FT(1) - FT(1) / FT(std::pow(10.0, 300.0)), FT(1) / FT(std::pow(10.0, 200.8))) ,
    Point_2(FT(1) - FT(1) / FT(std::pow(10.0, 300.0)), FT(1) / FT(std::pow(10.0, 200.5))) ,
    Point_2(FT(1) - FT(1) / FT(std::pow(10.0, 300.0)), FT(1) / FT(std::pow(10.0, 200.2))) ,
    Point_2(FT(1) - FT(1) / FT(std::pow(10.0, 300.0)), FT(1) / FT(std::pow(10.0, 200.05))),

    Point_2(
      FT(1) / FT(2),
      FT(1) / FT(2) * FT(vertices[0].y() + vertices[3].y()))
    };

  std::size_t count = 0;
  const Point_2 zero(0, 0);
  std::vector<FT> coordinates;
  for (std::size_t i = 0; i < 31; ++i) {

    CGAL::Barycentric_coordinates::discrete_harmonic_coordinates_2(
      vertices, queries[i], std::back_inserter(coordinates));
    const FT coordinate_sum =
      coordinates[count + 0] +
      coordinates[count + 1] +
      coordinates[count + 2] +
      coordinates[count + 3] ;
    const Point_2 linear_combination(
      vertices[0].x() * coordinates[count + 0] +
      vertices[1].x() * coordinates[count + 1] +
      vertices[2].x() * coordinates[count + 2] +
      vertices[3].x() * coordinates[count + 3] ,
      vertices[0].y() * coordinates[count + 0] +
      vertices[1].y() * coordinates[count + 1] +
      vertices[2].y() * coordinates[count + 2] +
      vertices[3].y() * coordinates[count + 3] );
    const Point_2 difference(
      linear_combination.x() - queries[i].x(),
      linear_combination.y() - queries[i].y());
    assert( (coordinate_sum == FT(1)) && (difference == zero) );
    count += 4;
  }

  const FT quater = FT(1) / FT(4);
  assert(
    coordinates[120] - quater == FT(0) &&
    coordinates[121] - quater == FT(0) &&
    coordinates[122] - quater == FT(0) &&
    coordinates[123] - quater == FT(0) );

  std::cout << "test_dh_almost_degenerate_polygon: PASSED" << std::endl;
  return EXIT_SUCCESS;
}
