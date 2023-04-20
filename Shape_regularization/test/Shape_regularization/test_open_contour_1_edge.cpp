#include "include/Saver.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Shape_regularization/regularize_contours.h>

namespace SR = CGAL::Shape_regularization;

template<class Traits>
void test_open_contour_1_edge() {

  using Point_2 = typename Traits::Point_2;
  using Contour = std::vector<Point_2>;
  using Saver   = SR::Tests::Saver<Traits>;

  using CD = SR::Contours::Longest_direction_2<Traits, Contour>;

  Saver saver;
  const Contour contour = {
    Point_2(0, 0), Point_2(1, 0)
  };
  assert(contour.size() == 2);

  // saver.export_open_contour(contour, "op1_input", 100);

  const bool is_closed = false;
  CD directions(contour, is_closed);

  std::vector<Point_2> regularized;
  SR::Contours::regularize_open_contour(
    contour, directions, std::back_inserter(regularized));

  const std::size_t num_directions =
    directions.number_of_directions();

  // saver.export_open_contour(regularized, "op1_output", 100);

  assert(num_directions == 1);
  assert(regularized.size() == 0);
}

int main() {
  test_open_contour_1_edge< CGAL::Simple_cartesian<double> >();
  test_open_contour_1_edge< CGAL::Exact_predicates_inexact_constructions_kernel >();
  test_open_contour_1_edge< CGAL::Exact_predicates_exact_constructions_kernel >();
  std::cout << "test_open_contour_1_edge: SUCCESS" << std::endl;
  return EXIT_SUCCESS;
}
