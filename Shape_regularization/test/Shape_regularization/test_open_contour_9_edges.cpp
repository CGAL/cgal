#include "include/Saver.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Shape_regularization/regularize_contours.h>

namespace SR = CGAL::Shape_regularization;

template<class Traits>
void test_open_contour_9_edges() {

  using Point_2 = typename Traits::Point_2;
  using Contour = std::vector<Point_2>;
  using Saver   = SR::Tests::Saver<Traits>;

  using CD = SR::Contours::Longest_direction_2<Traits, Contour>;

  Saver saver;
  const Contour contour = {
    Point_2(0.0, 0.0),
    Point_2(4.0, 0.0),
    Point_2(3.815571, 1.503828),
    Point_2(4.518233, 1.605529),
    Point_2(4.0, 2.5),
    Point_2(4.305586, 2.992361),
    Point_2(4.305586, 3.990881),
    Point_2(2.0, 3.5),
    Point_2(0.0, 4.0),
    Point_2(0.182071, 0.505309),
  };
  assert(contour.size() == 10);

  // saver.export_open_contour(contour, "op9_input", 100);

  const bool is_closed = false;
  CD directions(contour, is_closed);

  std::vector<Point_2> regularized;
  SR::Contours::regularize_open_contour(
    contour, directions, std::back_inserter(regularized),
    CGAL::parameters::default_values());

  const std::size_t num_directions =
    directions.number_of_directions();

  // saver.export_open_contour(regularized, "op9_output", 100);

  assert(num_directions == 1);
  assert(regularized.size() == 7);
}

int main() {
  test_open_contour_9_edges< CGAL::Simple_cartesian<double> >();
  test_open_contour_9_edges< CGAL::Exact_predicates_inexact_constructions_kernel >();
  std::cout << "test_open_contour_9_edges: SUCCESS" << std::endl;
  return EXIT_SUCCESS;
}
