#include "include/Saver.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Shape_regularization/regularize_segments.h>

namespace SR = CGAL::Shape_regularization;

template<class Traits>
void test_groups_orthogonal() {

  using FT        = typename Traits::FT;
  using Point_2   = typename Traits::Point_2;
  using Segment_2 = typename Traits::Segment_2;
  using Segments  = std::vector<Segment_2>;
  using Indices   = std::vector<std::size_t>;
  using Saver     = SR::Tests::Saver<Traits>;

  Saver saver;
  const Segments segments = {
    Segment_2(Point_2(2, 4), Point_2(2, 1)), // the left group
    Segment_2(Point_2(4, 1), Point_2(FT(39) / FT(10), 4)),
    Segment_2(Point_2(2, 6), Point_2(7, 6)),
    Segment_2(Point_2(5, 5), Point_2(5, 4)),
    Segment_2(Point_2(5, 3), Point_2(7, 3)),

    Segment_2(Point_2(7, 1), Point_2(10, 0)), // the bottom right group

    Segment_2(Point_2( 8, 3), Point_2(10, FT(52) / FT(10))), // the top right group
    Segment_2(Point_2(12, 2), Point_2(10, 4))
  };

  // saver.export_eps_segments(segments, "go_input", 100);

  std::vector<Indices> groups;
  SR::Segments::orthogonal_groups(
    segments, std::back_inserter(groups));
  assert(groups.size() == 3);

  // saver.export_eps_group(segments, groups[0], "go_group0", 100);
  assert(groups[0].size() == 5);
  assert(groups[0][0] == 0);
  assert(groups[0][1] == 1);
  assert(groups[0][2] == 3);
  assert(groups[0][3] == 2);
  assert(groups[0][4] == 4);
  // saver.export_eps_group(segments, groups[1], "go_group1", 100);
  assert(groups[1].size() == 1);
  assert(groups[1][0] == 5);
  // saver.export_eps_group(segments, groups[2], "go_group2", 100);
  assert(groups[2].size() == 2);
  assert(groups[2][0] == 6);
  assert(groups[2][1] == 7);

  groups.clear();
    SR::Segments::orthogonal_groups(
    segments, std::back_inserter(groups), CGAL::parameters::preserve_order(true));
  assert(groups.size() == 3);

  // saver.export_eps_group(segments, groups[0], "go_group0_pr", 100);
  assert(groups[0].size() == 5);
  assert(groups[0][0] == 0);
  assert(groups[0][1] == 1);
  assert(groups[0][2] == 2);
  assert(groups[0][3] == 3);
  assert(groups[0][4] == 4);
  // saver.export_eps_group(segments, groups[1], "go_group1_pr", 100);
  assert(groups[1].size() == 1);
  assert(groups[1][0] == 5);
  // saver.export_eps_group(segments, groups[2], "go_group2_pr", 100);
  assert(groups[2].size() == 2);
  assert(groups[2][0] == 6);
  assert(groups[2][1] == 7);
}

int main() {
  test_groups_orthogonal< CGAL::Simple_cartesian<double> >();
  test_groups_orthogonal< CGAL::Exact_predicates_inexact_constructions_kernel >();
  test_groups_orthogonal< CGAL::Exact_predicates_exact_constructions_kernel >();
  std::cout << "test_groups_orthogonal: SUCCESS" << std::endl;
  return EXIT_SUCCESS;
}
