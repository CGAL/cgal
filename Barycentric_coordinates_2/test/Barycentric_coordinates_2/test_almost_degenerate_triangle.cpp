#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/triangle_coordinates_2.h>

using Kernel  = CGAL::Exact_predicates_exact_constructions_kernel;
using FT      = typename Kernel::FT;
using Point_2 = typename Kernel::Point_2;

int main() {

  const Point_2 p0 = Point_2(0, 0);
  const Point_2 p1 = Point_2(1, 0);
  const Point_2 p2 = Point_2(-FT(1) / FT(2), FT(1) / FT(std::pow(10.0, 200.0)));

  const Point_2 queries[12] = {
    Point_2(0, FT(1) / FT(std::pow(10.0, 201.0))),
    Point_2(0, FT(1) / FT(std::pow(10.0, 220.0))),
    Point_2(0, FT(1) / FT(std::pow(10.0, 240.0))),
    Point_2(0, FT(1) / FT(std::pow(10.0, 260.0))),
    Point_2(0, FT(1) / FT(std::pow(10.0, 280.0))),
    Point_2(0, FT(1) / FT(std::pow(10.0, 300.0))),
    Point_2(FT(1) / FT(2), FT(1) / FT(std::pow(10.0, 200.5))),
    Point_2(FT(1) / FT(2), FT(1) / FT(std::pow(10.0, 201.0))),
    Point_2(FT(1) / FT(2), FT(1) / FT(std::pow(10.0, 260.0))),
    Point_2(FT(1) / FT(2), FT(1) / FT(std::pow(10.0, 280.0))),
    Point_2(FT(1) / FT(2), FT(1) / FT(std::pow(10.0, 300.0))),
    Point_2(
      (FT(1) / FT(3)) * FT(p0.x() + p1.x() + p2.x()) ,
      (FT(1) / FT(3)) * FT(p0.y() + p1.y() + p2.y()) )
    };

  std::size_t count = 0;
  const Point_2 zero(0, 0);
  std::vector<FT> coordinates;
  for (std::size_t i = 0; i < 12; ++i) {

    CGAL::Barycentric_coordinates::triangle_coordinates_2(
      p0, p1, p2, queries[i], std::back_inserter(coordinates));
    Point_2 linear_combination = Point_2(
      p0.x() * coordinates[count + 0] +
      p1.x() * coordinates[count + 1] +
      p2.x() * coordinates[count + 2] ,
      p0.y() * coordinates[count + 0] +
      p1.y() * coordinates[count + 1] +
      p2.y() * coordinates[count + 2] );
    Point_2 difference = Point_2(
      linear_combination.x() - queries[i].x(),
      linear_combination.y() - queries[i].y());
    assert(difference == zero);
    count += 3;

    const auto tuple =
      CGAL::Barycentric_coordinates::triangle_coordinates_in_tuple_2(
        p0, p1, p2, queries[i]);
    linear_combination = Point_2(
      p0.x() * std::get<0>(tuple) +
      p1.x() * std::get<1>(tuple) +
      p2.x() * std::get<2>(tuple) ,
      p0.y() * std::get<0>(tuple) +
      p1.y() * std::get<1>(tuple) +
      p2.y() * std::get<2>(tuple) );
    difference = Point_2(
      linear_combination.x() - queries[i].x(),
      linear_combination.y() - queries[i].y());
    assert(difference == zero);
  }

  const FT third = FT(1) / FT(3);
  assert(
    coordinates[33] - third == FT(0) &&
    coordinates[34] - third == FT(0) &&
    coordinates[35] - third == FT(0));

  std::cout << "test_almost_degenerate_triangle: PASSED" << std::endl;
  return EXIT_SUCCESS;
}
