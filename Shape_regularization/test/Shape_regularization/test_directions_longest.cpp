#include "include/Saver.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Shape_regularization/Contours/Longest_direction_2.h>

namespace SR = CGAL::Shape_regularization;

template<class Traits>
void test_directions_longest() {

  using Point_2     = typename Traits::Point_2;
  using Direction_2 = typename Traits::Direction_2;
  using Contour     = std::vector<Point_2>;
  using Saver       = SR::Tests::Saver<Traits>;

  using LD = SR::Contours::Longest_direction_2<Traits, Contour>;

  Saver saver;;
  const Contour contour = {
    Point_2(0, 0), Point_2(2, 0),
    Point_2(2, 1), Point_2(1, 1)
  };
  assert(contour.size() == 4);

  // saver.export_closed_contour(contour, "dl_input", 100);

  const bool is_closed = true;
  LD closed_directions(contour,  is_closed);
  LD   open_directions(contour, !is_closed);

  const std::size_t num_closed_directions =
    closed_directions.number_of_directions();
  const std::size_t num_open_directions =
    open_directions.number_of_directions();

  assert(num_closed_directions == 1);
  assert(num_closed_directions == num_open_directions);

  const auto& closed_dirs = closed_directions.get_directions();
  const auto& open_dirs = open_directions.get_directions();

  assert(closed_dirs[0] == open_dirs[0]);
  assert(closed_dirs[0] == Direction_2(contour[1] - contour[0]));
}

int main() {
  test_directions_longest< CGAL::Simple_cartesian<double> >();
  test_directions_longest< CGAL::Exact_predicates_inexact_constructions_kernel >();
  test_directions_longest< CGAL::Exact_predicates_exact_constructions_kernel >();
  std::cout << "test_directions_longest: SUCCESS" << std::endl;
  return EXIT_SUCCESS;
}
