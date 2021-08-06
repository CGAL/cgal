#include "include/utils.h"
#include "include/Saver.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Shape_regularization/regularize_segments.h>

namespace SR = CGAL::Shape_regularization;

template<class Traits>
void test_100_segments_offsets() {

  using FT        = typename Traits::FT;
  using Segment_2 = typename Traits::Segment_2;
  using Segments  = std::vector<Segment_2>;
  using Indices   = std::vector<std::size_t>;
  using Saver     = SR::Tests::Saver<Traits>;

  using NQ = SR::Segments::Delaunay_neighbor_query_2<Traits, Segments>;
  using OR = SR::Segments::Offset_regularization_2<Traits, Segments>;

  Saver saver;
  Segments segments;
  SR::Tests::create_example_offsets(segments);
  assert(segments.size() == 100);

  // saver.export_segments(segments, "100o_input", 100);

  const FT max_angle_2 = FT(1);
  std::vector<Indices> parallel_groups;
  SR::Segments::parallel_groups(
    segments,
    std::back_inserter(parallel_groups),
    CGAL::parameters::maximum_angle(max_angle_2));

  // Segments output;
  // for (std::size_t i = 0; i < parallel_groups.size(); ++i) {
  //   const auto& parallel_group = parallel_groups[i];
  //   output.clear();
  //   for (const std::size_t idx : parallel_group)
  //     output.push_back(segments[idx]);
  //   saver.export_segments(output, "output_" + std::to_string(i), 100);
  // }

  assert(segments.size() == 100);
  assert(parallel_groups.size() == 25);

  const FT max_offset_2 = FT(1) / FT(4);
  OR offset_regularization(
    segments, CGAL::parameters::maximum_offset(max_offset_2));

  NQ neighbor_query(segments);
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

  const std::size_t num_segments_offsets =
    offset_regularization.number_of_modified_segments();

  // saver.export_segments(segments, "100o_offsets", 100);
  // saver.export_segments(unique_segments, "100o_unique", 100);

  assert(segments.size() == 100);
  assert(collinear_groups.size() == parallel_groups.size());
  assert(unique_segments.size() == collinear_groups.size());
  assert(num_segments_offsets == 100);
}

int main() {
  test_100_segments_offsets< CGAL::Simple_cartesian<double> >();
  test_100_segments_offsets< CGAL::Exact_predicates_inexact_constructions_kernel >();
  test_100_segments_offsets< CGAL::Exact_predicates_exact_constructions_kernel >();
  std::cout << "test_100_segments_offsets: SUCCESS" << std::endl;
  return EXIT_SUCCESS;
}
