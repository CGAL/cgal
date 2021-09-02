#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/segment_coordinates_2.h>

using Kernel  = CGAL::Exact_predicates_exact_constructions_kernel;
using FT      = typename Kernel::FT;
using Point_2 = typename Kernel::Point_2;

int main() {

  const Point_2 p0 = Point_2(0, 0);
  const Point_2 p1 = Point_2(FT(1) / FT(std::pow(10.0, 200.0)), 0);

  const Point_2 queries[7] = {
    Point_2(FT(1) / FT(std::pow(10.0, 200.02)), 0),
    Point_2(FT(1) / FT(std::pow(10.0, 200.05)), 0),
    Point_2(FT(1) / FT(std::pow(10.0, 201.0)),  0),
    Point_2(FT(1) / FT(std::pow(10.0, 220.0)),  0),
    Point_2(FT(1) / FT(std::pow(10.0, 230.0)),  0),
    Point_2(FT(1) / FT(std::pow(10.0, 260.0)),  0),
    Point_2(
      (FT(1) / FT(2)) * FT(p0.x() + p1.x()) ,
      (FT(1) / FT(2)) * FT(p0.y() + p1.y()) )
  };

  std::size_t count = 0;
  const Point_2 zero(0, 0);
  std::vector<FT> coordinates;
  for (std::size_t i = 0; i < 7; ++i) {

    CGAL::Barycentric_coordinates::segment_coordinates_2(
      p0, p1, queries[i], std::back_inserter(coordinates));
    Point_2 linear_combination = Point_2(
      p0.x() * coordinates[count + 0] + p1.x() * coordinates[count + 1],
      p0.y() * coordinates[count + 0] + p1.y() * coordinates[count + 1]);
    Point_2 difference = Point_2(
      linear_combination.x() - queries[i].x(),
      linear_combination.y() - queries[i].y());
    assert(difference == zero);
    count += 2;

    const auto pair =
      CGAL::Barycentric_coordinates::segment_coordinates_in_pair_2(
        p0, p1, queries[i]);
    linear_combination = Point_2(
      p0.x() * pair.first + p1.x() * pair.second,
      p0.y() * pair.first + p1.y() * pair.second);
    difference = Point_2(
      linear_combination.x() - queries[i].x(),
      linear_combination.y() - queries[i].y());
    assert(difference == zero);
  }

  const FT half = FT(1) / FT(2);
  assert(
    coordinates[12] - half == FT(0) &&
    coordinates[13] - half == FT(0));

  std::cout << "test_almost_degenerate_segment: PASSED" << std::endl;
  return EXIT_SUCCESS;
}
