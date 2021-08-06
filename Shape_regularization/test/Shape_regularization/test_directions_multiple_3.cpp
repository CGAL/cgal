#include "include/utils.h"
#include "include/Saver.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Shape_regularization/Contours/Multiple_directions_2.h>
#include <CGAL/Shape_regularization/regularize_contours.h>

namespace SR = CGAL::Shape_regularization;

template<class Traits>
void test_directions_multiple_3() {

  using Point_2 = typename Traits::Point_2;
  using Contour = std::vector<Point_2>;
  using Saver   = SR::Tests::Saver<Traits>;

  using MD = SR::Contours::Multiple_directions_2<Traits, Contour>;

  Saver saver;
  Contour contour;
  SR::Tests::create_example_multiple_directions(contour);
  assert(contour.size() == 29);

  // saver.export_closed_contour(contour, "dm3_input", 100);

  const bool is_closed = true;
  MD closed_directions(
    contour,  is_closed, CGAL::parameters::adjust_directions(false));
  MD open_directions(
    contour, !is_closed, CGAL::parameters::adjust_directions(false));

  const std::size_t num_closed_directions =
    closed_directions.number_of_directions();
  const std::size_t num_open_directions =
    open_directions.number_of_directions();

  assert(num_closed_directions == 3);
  assert(num_closed_directions == num_open_directions);

  const auto& closed_dirs = closed_directions.get_directions();
  const auto& open_dirs = open_directions.get_directions();

  // std::vector<Point_2> regularized;
  // SR::Contours::regularize_closed_contour(
  //   contour, closed_directions, std::back_inserter(regularized),
  //   CGAL::parameters::maximum_offset(2));
  // saver.export_closed_contour(regularized, "dm3_output_cl", 100);

  // regularized.clear();
  // SR::Contours::regularize_open_contour(
  //   contour, open_directions, std::back_inserter(regularized),
  //   CGAL::parameters::maximum_offset(2));
  // saver.export_open_contour(regularized, "dm3_output_op", 100);

  assert(closed_dirs[0] == open_dirs[0]);
  assert(closed_dirs[1] == open_dirs[1]);
  assert(closed_dirs[2] == open_dirs[2]);
}

int main() {
  test_directions_multiple_3< CGAL::Simple_cartesian<double> >();
  test_directions_multiple_3< CGAL::Exact_predicates_inexact_constructions_kernel >();
  test_directions_multiple_3< CGAL::Exact_predicates_exact_constructions_kernel >();
  std::cout << "test_directions_multiple_3: SUCCESS" << std::endl;
  return EXIT_SUCCESS;
}
