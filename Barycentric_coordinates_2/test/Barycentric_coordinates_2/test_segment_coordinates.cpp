#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/segment_coordinates_2.h>

using Kernel  = CGAL::Exact_predicates_exact_constructions_kernel;
using FT      = typename Kernel::FT;
using Point_2 = typename Kernel::Point_2;

int main() {

  const Point_2 p0 = Point_2(0, 0);
  const Point_2 p1 = Point_2(1, 0);

  const Point_2 queries[6] = {
    Point_2( FT(2) / FT(5),  0),
    Point_2(1, 0),
    Point_2( FT(7) / FT(10), 0),
    Point_2(-FT(3) / FT(10), 0),
    Point_2( FT(6) / FT(5),  0),
    Point_2(0, 0)
  };

  const FT expected[12] = {
    FT(3)  / FT(5) ,  FT(2) / FT(5) ,
    0, 1,
    FT(3)  / FT(10),  FT(7) / FT(10),
    FT(13) / FT(10), -FT(3) / FT(10),
   -FT(1)  / FT(5) ,  FT(6) / FT(5) ,
    1, 0
  };

  std::size_t count = 0;
  std::vector<FT> coordinates;
  for (std::size_t i = 0; i < 6; ++i) {
    CGAL::Barycentric_coordinates::segment_coordinates_2(
      p0, p1, queries[i], std::back_inserter(coordinates));
    assert(
      coordinates[count + 0] - expected[count + 0] == FT(0) &&
      coordinates[count + 1] - expected[count + 1] == FT(0) );
    count += 2;
  }

  count = 0;
  for (std::size_t i = 0; i < 6; ++i) {
    const auto pair =
      CGAL::Barycentric_coordinates::segment_coordinates_in_pair_2(
        p0, p1, queries[i]);
    assert(
      pair.first  - expected[count + 0] == FT(0) &&
      pair.second - expected[count + 1] == FT(0) );
    count += 2;
  }

  std::cout << "test_segment_coordinates: PASSED" << std::endl;
  return EXIT_SUCCESS;
}
