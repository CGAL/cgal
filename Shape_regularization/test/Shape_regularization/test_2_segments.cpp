#include "include/utils.h"
#include "include/Saver.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Shape_regularization.h>

namespace SR = CGAL::Shape_regularization;

template<class Traits>
void test_2_segments() {

  using FT        = typename Traits::FT;
  using Point_2   = typename Traits::Point_2;
  using Segment_2 = typename Traits::Segment_2;
  using Segments  = std::vector<Segment_2>;
  using Indices   = std::vector<std::size_t>;
  using Saver     = SR::Tests::Saver<Traits>;

  using Segment_map = CGAL::Identity_property_map<Segment_2>;

  using NQ = SR::Segments::Delaunay_neighbor_query_2<Traits, Segments, Segment_map>;
  using AR = SR::Segments::Angle_regularization_2<Traits, Segments, Segment_map>;
  using OR = SR::Segments::Offset_regularization_2<Traits, Segments, Segment_map>;
  using QP = CGAL::OSQP_quadratic_program_traits<FT>;

  using QP_AR = SR::QP_regularization<Traits, Segments, NQ, AR, QP>;
  using QP_OR = SR::QP_regularization<Traits, Segments, NQ, OR, QP>;

  Saver saver;
  Segment_map smap;
  Segments segments = {
    Segment_2(Point_2(1            , 1), Point_2(1            , 4)),
    Segment_2(Point_2(FT(3) / FT(2), 4), Point_2(FT(3) / FT(2), 5))
  };
  assert(segments.size() == 2);

  // saver.export_segments(segments, "2_input", 100);

  const FT max_angle_2 = FT(5);
  NQ neighbor_query(
    segments, CGAL::parameters::
    segment_map(smap));
  AR angle_regularization(
    segments, CGAL::parameters::
    maximum_angle(max_angle_2).
    segment_map(smap));

  QP qp_angles;
  QP_AR qp_ar(
    segments, neighbor_query, angle_regularization, qp_angles);
  qp_ar.regularize();

  std::vector<Indices> parallel_groups;
  angle_regularization.parallel_groups(
    std::back_inserter(parallel_groups));

  std::vector<Indices> orthogonal_groups;
  angle_regularization.orthogonal_groups(
    std::back_inserter(orthogonal_groups));

  const std::size_t num_segments_angles =
    angle_regularization.number_of_modified_segments();

  // saver.export_segments(segments, "2_angles", 100);

  assert(segments.size() == 2);
  assert(parallel_groups.size() == 1);
  assert(orthogonal_groups.size() == 1);
  assert(num_segments_angles == 2);

  std::vector<int> reference_values;
  reference_values.reserve(2);
  reference_values.push_back(7);
  reference_values.push_back(12);
  SR::Tests::check_reference_values(segments, reference_values);

  const FT max_offset_2 = FT(1) / FT(2);
  OR offset_regularization(
    segments, CGAL::parameters::
    maximum_offset(max_offset_2).
    segment_map(smap));

  neighbor_query.clear();
  for (const auto& parallel_group : parallel_groups) {
    neighbor_query.add_group(parallel_group);
    offset_regularization.add_group(parallel_group);
  }

  QP qp_offsets;
  QP_OR qp_or(
    segments, neighbor_query, offset_regularization, qp_offsets);
  qp_or.regularize();

  std::vector<Indices> collinear_groups;
  offset_regularization.collinear_groups(
    std::back_inserter(collinear_groups));

  std::vector<Segment_2> unique_segments;
  offset_regularization.unique_segments(
    std::back_inserter(unique_segments));

  const std::size_t num_segments_offsets =
    offset_regularization.number_of_modified_segments();

  // saver.export_segments(segments, "2_offsets", 100);
  // saver.export_segments(unique_segments, "2_unique", 100);

  assert(segments.size() == 2);
  assert(collinear_groups.size() == 1);
  assert(collinear_groups[0].size() == 2);
  assert(unique_segments.size() == 1);
  assert(num_segments_offsets == 2);

  reference_values.clear();
  reference_values.push_back(7);
  reference_values.push_back(11);
  SR::Tests::check_reference_values(segments, reference_values);
}

int main() {
  test_2_segments< CGAL::Simple_cartesian<double> >();
  test_2_segments< CGAL::Exact_predicates_inexact_constructions_kernel >();
  test_2_segments< CGAL::Exact_predicates_exact_constructions_kernel >();
  std::cout << "test_2_segments: SUCCESS" << std::endl;
  return EXIT_SUCCESS;
}
