#include "include/Saver.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Shape_regularization/regularize_segments.h>

namespace SR = CGAL::Shape_regularization;

template<class Traits>
void test_groups_collinear() {

  using FT        = typename Traits::FT;
  using Point_2   = typename Traits::Point_2;
  using Segment_2 = typename Traits::Segment_2;
  using Segments  = std::vector<Segment_2>;
  using Indices   = std::vector<std::size_t>;
  using Saver     = SR::Tests::Saver<Traits>;

  Saver saver;
  const Segments segments = {
    Segment_2(Point_2(2, -1), Point_2(5, -1)), // the bottom group
    Segment_2(Point_2(6, -1), Point_2(7, -1)),

    Segment_2(Point_2(1, 1), Point_2(4, 1)), // the middle group
    Segment_2(Point_2(5, FT(11) / FT(10)), Point_2(6, FT(11) / FT(10))),
    Segment_2(Point_2(4, FT(9)  / FT(10)), Point_2(6, FT(9)  / FT(10))),

    Segment_2(Point_2(2, 2), Point_2(3, 3)), // the top left group
    Segment_2(Point_2(5, 5), Point_2(4, 4)),

    Segment_2(Point_2(7, 2), Point_2(5, 4)) // the top right group
  };

  // saver.export_eps_segments(segments, "gc_input", 100);

  std::vector<Indices> groups;
  SR::Segments::collinear_groups(
    segments, std::back_inserter(groups));
  assert(groups.size() == 4);

  // saver.export_eps_group(segments, groups[0], "gc_group0", 100);
  assert(groups[0].size() == 2);
  assert(groups[0][0] == 0);
  assert(groups[0][1] == 1);
  // saver.export_eps_group(segments, groups[1], "gc_group1", 100);
  assert(groups[1].size() == 3);
  assert(groups[1][0] == 2);
  assert(groups[1][1] == 3);
  assert(groups[1][2] == 4);
  // saver.export_eps_group(segments, groups[2], "gc_group2", 100);
  assert(groups[2].size() == 2);
  assert(groups[2][0] == 5);
  assert(groups[2][1] == 6);
  // saver.export_eps_group(segments, groups[3], "gc_group3", 100);
  assert(groups[3].size() == 1);
  assert(groups[3][0] == 7);

  groups.clear();
  SR::Segments::collinear_groups(
    segments, std::back_inserter(groups), CGAL::parameters::preserve_order(true));
  assert(groups.size() == 4);

  // saver.export_eps_group(segments, groups[0], "gc_group0_pr", 100);
  assert(groups[0].size() == 2);
  assert(groups[0][0] == 0);
  assert(groups[0][1] == 1);
  // saver.export_eps_group(segments, groups[1], "gc_group1_pr", 100);
  assert(groups[1].size() == 3);
  assert(groups[1][0] == 2);
  assert(groups[1][1] == 3);
  assert(groups[1][2] == 4);
  // saver.export_eps_group(segments, groups[2], "gc_group2_pr", 100);
  assert(groups[2].size() == 2);
  assert(groups[2][0] == 5);
  assert(groups[2][1] == 6);
  // saver.export_eps_group(segments, groups[3], "gc_group3_pr", 100);
  assert(groups[3].size() == 1);
  assert(groups[3][0] == 7);
}

int main() {
  test_groups_collinear< CGAL::Simple_cartesian<double> >();
  test_groups_collinear< CGAL::Exact_predicates_inexact_constructions_kernel >();
  test_groups_collinear< CGAL::Exact_predicates_exact_constructions_kernel >();
  std::cout << "test_groups_collinear: SUCCESS" << std::endl;
  return EXIT_SUCCESS;
}
