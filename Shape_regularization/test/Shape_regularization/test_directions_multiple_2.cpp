#include "include/Saver.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Shape_regularization/Contours/Multiple_directions_2.h>
#include <CGAL/Shape_regularization/regularize_contours.h>

namespace SR = CGAL::Shape_regularization;

template<class Traits>
void test_directions_multiple_2() {

  using FT          = typename Traits::FT;
  using Point_2     = typename Traits::Point_2;
  using Direction_2 = typename Traits::Direction_2;
  using Contour     = std::vector<Point_2>;
  using Saver       = SR::Tests::Saver<Traits>;

  using Point_map = CGAL::Identity_property_map<Point_2>;
  using MD = SR::Contours::Multiple_directions_2<Traits, Contour, Point_map>;

  Saver saver;
  Point_map pmap;
  const Contour contour = {
    Point_2( 1, 1), Point_2(4, 1),
    Point_2( 4, 4), Point_2(7, 1),
    Point_2(10, 4), Point_2(7, 7),
    Point_2(1, 7)
  };
  assert(contour.size() == 7);

  // saver.export_closed_contour(contour, "dm2_input", 100);

  const FT min_length_2 = FT(2);
  const FT max_angle_2 = FT(10);

  const bool is_closed = true;
  MD closed_directions(
    contour,  is_closed,
    CGAL::parameters::
    minimum_length(min_length_2).
    maximum_angle(max_angle_2).
    adjust_directions(false).
    point_map(pmap));
  MD open_directions(
    contour, !is_closed,
    CGAL::parameters::
    minimum_length(min_length_2).
    maximum_angle(max_angle_2).
    adjust_directions(false).
    point_map(pmap));

  const std::size_t num_closed_directions =
    closed_directions.number_of_directions();
  const std::size_t num_open_directions =
    open_directions.number_of_directions();

  assert(num_closed_directions == 2);
  assert(num_closed_directions == num_open_directions);

  const auto& closed_dirs = closed_directions.get_directions();
  const auto& open_dirs = open_directions.get_directions();

  const auto d = contour[5] - contour[4];
  const auto length = static_cast<FT>(
    CGAL::sqrt(CGAL::to_double(d.x() * d.x() + d.y() * d.y())));
  const std::vector<Direction_2> ref_directions = {
    Direction_2(1, 0),
    Direction_2(d.x() / length, d.y() / length)
  };

  // std::vector<Point_2> regularized;
  // SR::Contours::regularize_closed_contour(
  //   contour, closed_directions, std::back_inserter(regularized),
  //   CGAL::parameters::maximum_offset(2));
  // saver.export_closed_contour(regularized, "dm2_output_cl", 100);

  // regularized.clear();
  // SR::Contours::regularize_open_contour(
  //   contour, open_directions, std::back_inserter(regularized),
  //   CGAL::parameters::maximum_offset(2));
  // saver.export_open_contour(regularized, "dm2_output_op", 100);

  assert(closed_dirs[0] == open_dirs[0]);
  assert(closed_dirs[1] == open_dirs[1]);
  assert(closed_dirs[0] == ref_directions[0]);
  assert(closed_dirs[1] == ref_directions[1]);
}

int main() {
  test_directions_multiple_2< CGAL::Simple_cartesian<double> >();
  test_directions_multiple_2< CGAL::Exact_predicates_inexact_constructions_kernel >();
  test_directions_multiple_2< CGAL::Exact_predicates_exact_constructions_kernel >();
  std::cout << "test_directions_multiple_2: SUCCESS" << std::endl;
  return EXIT_SUCCESS;
}
