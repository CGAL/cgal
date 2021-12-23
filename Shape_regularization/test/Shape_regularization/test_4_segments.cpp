#include "include/utils.h"
#include "include/Saver.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Shape_regularization/regularize_segments.h>

namespace SR = CGAL::Shape_regularization;

template<class Traits>
void test_4_segments() {

  using FT        = typename Traits::FT;
  using Point_2   = typename Traits::Point_2;
  using Segment_2 = typename Traits::Segment_2;
  using Segments  = std::vector<Segment_2>;
  using Indices   = std::vector<std::size_t>;
  using Saver     = SR::Tests::Saver<Traits>;

  using NQ = SR::Segments::Delaunay_neighbor_query_2<Traits, Segments>;
  using AR = SR::Segments::Angle_regularization_2<Traits, Segments>;
  using OR = SR::Segments::Offset_regularization_2<Traits, Segments>;
  using QP = CGAL::OSQP_quadratic_program_traits<FT>;

  Saver saver;
  Segments segments = {
    Segment_2(Point_2(FT(2) / FT(10), 1)              , Point_2(FT(2)  / FT(10) , FT(2)  / FT(10))),
    Segment_2(Point_2(FT(2) / FT(10), 1)              , Point_2(FT(25) / FT(100), FT(16) / FT(10))),
    Segment_2(Point_2(FT(6) / FT(10), FT(16) / FT(10)), Point_2(FT(6)  / FT(10) , FT(14) / FT(10))),
    Segment_2(Point_2(FT(2) / FT(10), 2)              , Point_2(FT(6)  / FT(10) , 2))
  };
  assert(segments.size() == 4);

  // saver.export_segments(segments, "4_input", 100);

  const FT max_angle_2 = FT(5);
  NQ neighbor_query(segments);
  AR angle_regularization(
    segments, CGAL::parameters::maximum_angle(max_angle_2));

  QP qp_angles;
  SR::Segments::regularize_segments(
    segments, neighbor_query, angle_regularization, qp_angles);

  std::vector<Indices> parallel_groups;
  angle_regularization.parallel_groups(
    std::back_inserter(parallel_groups));

  std::vector<Indices> orthogonal_groups;
  angle_regularization.orthogonal_groups(
    std::back_inserter(orthogonal_groups));

  const std::size_t num_segments_angles =
    angle_regularization.number_of_modified_segments();

  // saver.export_segments(segments, "4_angles", 100);

  assert(segments.size() == 4);
  assert(parallel_groups.size() == 2);
  assert(orthogonal_groups.size() == 1);
  assert(num_segments_angles == 4);

  std::vector<int> reference_values;
  reference_values.reserve(4);
  reference_values.push_back(1);
  reference_values.push_back(3);
  reference_values.push_back(4);
  reference_values.push_back(4);
  SR::Tests::check_reference_values(segments, reference_values);

  const FT max_offset_2 = FT(1) / FT(10);
  OR offset_regularization(
    segments, CGAL::parameters::maximum_offset(max_offset_2));

  neighbor_query.clear();
  for (const auto& parallel_group : parallel_groups) {
    neighbor_query.add_group(parallel_group);
    offset_regularization.add_group(parallel_group);
  }

  QP qp_offsets;
  SR::Segments::regularize_segments(
    segments, neighbor_query, offset_regularization, qp_offsets);

  std::vector<Indices> collinear_groups;
  offset_regularization.collinear_groups(
    std::back_inserter(collinear_groups));

  std::vector<Segment_2> unique_segments;
  offset_regularization.unique_segments(
    std::back_inserter(unique_segments));

  const std::size_t num_segments_offsets =
    offset_regularization.number_of_modified_segments();

  // saver.export_segments(segments, "4_offsets", 100);
  // saver.export_segments(unique_segments, "4_unique", 100);

  assert(segments.size() == 4);
  assert(collinear_groups.size() == 3);
  assert(collinear_groups[0].size() == 2);
  assert(collinear_groups[1].size() == 1);
  assert(collinear_groups[2].size() == 1);
  assert(unique_segments.size() == 3);
  assert(num_segments_offsets == 3);

  reference_values.clear();
  reference_values.push_back(1);
  reference_values.push_back(3);
  reference_values.push_back(4);
  reference_values.push_back(4);
  SR::Tests::check_reference_values(segments, reference_values);
}

int main() {
  test_4_segments< CGAL::Simple_cartesian<double> >();
  test_4_segments< CGAL::Exact_predicates_inexact_constructions_kernel >();
  test_4_segments< CGAL::Exact_predicates_exact_constructions_kernel >();
  std::cout << "test_4_segments: SUCCESS" << std::endl;
  return EXIT_SUCCESS;
}
