#include "include/Saver.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Shape_regularization/regularize_segments.h>

namespace SR = CGAL::Shape_regularization;

template<class Traits>
void test_groups_parallel() {

  using FT        = typename Traits::FT;
  using Point_2   = typename Traits::Point_2;
  using Segment_2 = typename Traits::Segment_2;
  using Segments  = std::vector<Segment_2>;
  using Indices   = std::vector<std::size_t>;
  using Saver     = SR::Tests::Saver<Traits>;

  Saver saver;
  const Segments segments = {
    Segment_2(Point_2(1, 1), Point_2(4, 1)), // the bottom group
    Segment_2(Point_2(1, 2), Point_2(4, 2)),
    Segment_2(Point_2(1, 3), Point_2(FT(399) / FT(100), FT(319) / FT(100))),

    Segment_2(Point_2(1, 4), Point_2(1, 6)), // the top left group
    Segment_2(Point_2(2, 5), Point_2(2, 8)),

    Segment_2(Point_2(3, 5), Point_2(6, 7)), // the top right group
    Segment_2(Point_2(7, 6), Point_2(4, 4))
  };

  // saver.export_eps_segments(segments, "gp_input", 100);

  std::vector<Indices> groups;
  SR::Segments::parallel_groups(
    segments, std::back_inserter(groups));
  assert(groups.size() == 3);

  // saver.export_eps_group(segments, groups[0], "gp_group0", 100);
  assert(groups[0].size() == 3);
  assert(groups[0][0] == 0);
  assert(groups[0][1] == 1);
  assert(groups[0][2] == 2);
  // saver.export_eps_group(segments, groups[1], "gp_group1", 100);
  assert(groups[1].size() == 2);
  assert(groups[1][0] == 3);
  assert(groups[1][1] == 4);
  // saver.export_eps_group(segments, groups[2], "gp_group2", 100);
  assert(groups[2].size() == 2);
  assert(groups[2][0] == 5);
  assert(groups[2][1] == 6);

  groups.clear();
  SR::Segments::parallel_groups(
    segments, std::back_inserter(groups), CGAL::parameters::preserve_order(true));
  assert(groups.size() == 3);

  // saver.export_eps_group(segments, groups[0], "gp_group0_pr", 100);
  assert(groups[0].size() == 3);
  assert(groups[0][0] == 0);
  assert(groups[0][1] == 1);
  assert(groups[0][2] == 2);
  // saver.export_eps_group(segments, groups[1], "gp_group1_pr", 100);
  assert(groups[1].size() == 2);
  assert(groups[1][0] == 3);
  assert(groups[1][1] == 4);
  // saver.export_eps_group(segments, groups[2], "gp_group2_pr", 100);
  assert(groups[2].size() == 2);
  assert(groups[2][0] == 5);
  assert(groups[2][1] == 6);
}

int main() {
  test_groups_parallel< CGAL::Simple_cartesian<double> >();
  test_groups_parallel< CGAL::Exact_predicates_inexact_constructions_kernel >();
  test_groups_parallel< CGAL::Exact_predicates_exact_constructions_kernel >();
  std::cout << "test_groups_parallel: SUCCESS" << std::endl;
  return EXIT_SUCCESS;
}
