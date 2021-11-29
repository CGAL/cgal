#include "include/Saver.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Shape_regularization/regularize_segments.h>

namespace SR = CGAL::Shape_regularization;

template<typename Segment_2>
bool has_no_difference(
  const Segment_2& a,
  const Segment_2& b) {

  const auto& as = a.source();
  const auto& at = a.target();
  const auto& bs = b.source();
  const auto& bt = b.target();

  const auto v1 = bs - as;
  const auto v2 = bt - at;
  const auto v3 = bs - at;
  const auto v4 = bt - as;

  const double l1 = CGAL::abs(CGAL::to_double(v1.x() * v1.x() + v1.y() * v1.y()));
  const double l2 = CGAL::abs(CGAL::to_double(v2.x() * v2.x() + v2.y() * v2.y()));
  const double l3 = CGAL::abs(CGAL::to_double(v3.x() * v3.x() + v3.y() * v3.y()));
  const double l4 = CGAL::abs(CGAL::to_double(v4.x() * v4.x() + v4.y() * v4.y()));

  const double tol = 0.000001;
  return (l1 < tol && l2 < tol) || (l3 < tol && l4 < tol);
}

template<class Traits>
void test_segments() {

  using Point_2   = typename Traits::Point_2;
  using Segment_2 = typename Traits::Segment_2;
  using Segments  = std::vector<Segment_2>;
  using Indices   = std::vector<std::size_t>;
  using Saver     = SR::Tests::Saver<Traits>;

  using NQ = SR::Segments::Delaunay_neighbor_query_2<Traits, Segments>;
  using AR = SR::Segments::Angle_regularization_2<Traits, Segments>;
  using OR = SR::Segments::Offset_regularization_2<Traits, Segments>;

  Saver saver;
  Segments segments = {
    Segment_2(Point_2(0, 0), Point_2(1, 0)), // bottom segments
    Segment_2(Point_2(0, 1), Point_2(1, 1)),

    Segment_2(Point_2(0, 2), Point_2(0, 3)), // top segments
    Segment_2(Point_2(0, 4), Point_2(0, 5)),
  };
  const auto saved = segments;
  assert(segments.size() == 4);

  // saver.export_segments(segments, "sg_input", 100);

  NQ neighbor_query(segments);
  AR angle_regularization(segments);
  SR::Segments::regularize_segments(
    segments, neighbor_query, angle_regularization);

  std::vector<Indices> parallel_groups;
  angle_regularization.parallel_groups(
    std::back_inserter(parallel_groups));

  std::vector<Indices> orthogonal_groups;
  angle_regularization.orthogonal_groups(
    std::back_inserter(orthogonal_groups));

  // saver.export_segments(segments, "sg_angles", 100);

  assert(segments.size() == 4);
  assert(has_no_difference(segments[0], saved[0]));
  assert(has_no_difference(segments[1], saved[1]));
  assert(has_no_difference(segments[2], saved[2]));
  assert(has_no_difference(segments[3], saved[3]));

  assert(parallel_groups.size() == 2);
  assert(parallel_groups[0].size() == 2);
  assert(parallel_groups[0][0] == 0);
  assert(parallel_groups[0][1] == 1);
  assert(parallel_groups[1].size() == 2);
  assert(parallel_groups[1][0] == 2);
  assert(parallel_groups[1][1] == 3);

  assert(orthogonal_groups.size() == 1);
  assert(orthogonal_groups[0].size() == 4);
  assert(orthogonal_groups[0][0] == 0);
  assert(orthogonal_groups[0][1] == 1);
  assert(orthogonal_groups[0][2] == 2);
  assert(orthogonal_groups[0][3] == 3);

  OR offset_regularization(segments);
  neighbor_query.clear();
  for (const auto& parallel_group : parallel_groups) {
    neighbor_query.add_group(parallel_group);
    offset_regularization.add_group(parallel_group);
  }

  SR::Segments::regularize_segments(
    segments, neighbor_query, offset_regularization);

  std::vector<Indices> collinear_groups;
  offset_regularization.collinear_groups(
    std::back_inserter(collinear_groups));

  std::vector<Segment_2> unique_segments;
  offset_regularization.unique_segments(
    std::back_inserter(unique_segments));

  // saver.export_segments(segments, "sg_offsets", 100);
  // saver.export_segments(unique_segments, "sg_unique", 100);

  assert(segments.size() == 4);
  assert(has_no_difference(segments[0], saved[0]));
  assert(has_no_difference(segments[1], saved[1]));
  assert(has_no_difference(segments[2], saved[2]));
  assert(has_no_difference(segments[3], saved[3]));

  assert(collinear_groups.size() == 3);
  assert(collinear_groups[0].size() == 1);
  assert(collinear_groups[0][0] == 0);
  assert(collinear_groups[1].size() == 1);
  assert(collinear_groups[1][0] == 1);
  assert(collinear_groups[2].size() == 2);
  assert(collinear_groups[2][0] == 2);
  assert(collinear_groups[2][1] == 3);

  assert(unique_segments.size() == 3);
  const Segment_2 ref_segment(Point_2(0, 2), Point_2(0, 5));
  assert(has_no_difference(unique_segments[0], saved[0]));
  assert(has_no_difference(unique_segments[1], saved[1]));
  assert(has_no_difference(unique_segments[2], ref_segment));
}

int main() {
  test_segments< CGAL::Simple_cartesian<double> >();
  test_segments< CGAL::Exact_predicates_inexact_constructions_kernel >();
  test_segments< CGAL::Exact_predicates_exact_constructions_kernel >();
  std::cout << "test_segments: SUCCESS" << std::endl;
  return EXIT_SUCCESS;
}
