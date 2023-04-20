#include "include/Saver.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Shape_regularization/regularize_contours.h>

namespace SR = CGAL::Shape_regularization;

template<class Traits>
void test_contour_orientation() {

  using Point_2 = typename Traits::Point_2;
  using Contour = std::vector<Point_2>;
  using Saver   = SR::Tests::Saver<Traits>;

  using CD = SR::Contours::Longest_direction_2<Traits, Contour>;

  Saver saver;
  const Contour direct = {
    Point_2(0.1, 0.0), Point_2(1.1, 0.0), Point_2(1.0, 1.0), Point_2(0.0, 1.0)
  };
  assert(direct.size() == 4);

  // saver.export_closed_contour(direct, "cl_direct_input", 100);
  // saver.export_open_contour(direct, "op_dir_input", 100);

  const Contour reversed = {
    Point_2(0.1, 0.0), Point_2(0.0, 1.0), Point_2(1.0, 1.0), Point_2(1.1, 0.0)
  };
  assert(reversed.size() == 4);

  // saver.export_closed_contour(reversed, "cl_reversed_input", 100);
  // saver.export_open_contour(reversed, "op_rev_input", 100);

  CD dclosed_direct(direct, true);
  CD dclosed_reversed(reversed, true);

  std::vector<Point_2> rclosed_direct, rclosed_reversed;
  SR::Contours::regularize_closed_contour(
    direct, dclosed_direct, std::back_inserter(rclosed_direct));
  SR::Contours::regularize_closed_contour(
    reversed, dclosed_reversed, std::back_inserter(rclosed_reversed));

  // saver.export_closed_contour(rclosed_direct, "cl_dir_output", 100);
  // saver.export_closed_contour(rclosed_reversed, "cl_rev_output", 100);

  assert(rclosed_direct.size() == 4);
  assert(rclosed_reversed.size() == 4);

  CD dopen_direct(direct, false);
  CD dopen_reversed(reversed, false);

  std::vector<Point_2> ropen_direct, ropen_reversed;
  SR::Contours::regularize_open_contour(
    direct, dopen_direct, std::back_inserter(ropen_direct));
  SR::Contours::regularize_open_contour(
    reversed, dopen_reversed, std::back_inserter(ropen_reversed));

  // saver.export_open_contour(ropen_direct, "op_dir_output", 100);
  // saver.export_open_contour(ropen_reversed, "op_rev_output", 100);

  assert(ropen_direct.size() == 4);
  assert(ropen_reversed.size() == 4);
}

int main() {
  test_contour_orientation< CGAL::Simple_cartesian<double> >();
  test_contour_orientation< CGAL::Exact_predicates_inexact_constructions_kernel >();
  test_contour_orientation< CGAL::Exact_predicates_exact_constructions_kernel >();
  std::cout << "test_contour_orientation: SUCCESS" << std::endl;
  return EXIT_SUCCESS;
}
