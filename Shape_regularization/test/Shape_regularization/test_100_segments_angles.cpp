#include "include/utils.h"
#include "include/Saver.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Shape_regularization/regularize_segments.h>

namespace SR = CGAL::Shape_regularization;

template<class Traits>
void test_100_segments_angles() {

  using FT        = typename Traits::FT;
  using Segment_2 = typename Traits::Segment_2;
  using Segments  = std::vector<Segment_2>;
  using Indices   = std::vector<std::size_t>;
  using Saver     = SR::Tests::Saver<Traits>;

  using NQ = SR::Segments::Delaunay_neighbor_query_2<Traits, Segments>;
  using AR = SR::Segments::Angle_regularization_2<Traits, Segments>;

  Saver saver;
  Segments segments;
  SR::Tests::create_example_angles(segments);
  assert(segments.size() == 100);

  // saver.export_segments(segments, "100a_input", 100);

  const FT max_angle_2 = FT(40);
  NQ neighbor_query(segments);
  AR angle_regularization(
    segments, CGAL::parameters::maximum_angle(max_angle_2));

  SR::Segments::regularize_segments(
    segments, neighbor_query, angle_regularization);

  std::vector<Indices> parallel_groups;
  angle_regularization.parallel_groups(
    std::back_inserter(parallel_groups));

  // Segments output;
  // for (std::size_t i = 0; i < parallel_groups.size(); ++i) {
  //   const auto& parallel_group = parallel_groups[i];
  //   output.clear();
  //   for (const std::size_t idx : parallel_group)
  //     output.push_back(segments[idx]);
  //   saver.export_segments(output, "output_" + std::to_string(i), 100);
  // }

  std::vector<Indices> orthogonal_groups;
  angle_regularization.orthogonal_groups(
    std::back_inserter(orthogonal_groups));

  const std::size_t num_segments_angles =
    angle_regularization.number_of_modified_segments();

  // saver.export_segments(segments, "100a_angles", 100);

  assert(segments.size() == 100);
  assert(parallel_groups.size() == 2);
  assert(orthogonal_groups.size() == 1);
  assert(num_segments_angles == 100);
}

int main() {
  test_100_segments_angles< CGAL::Simple_cartesian<double> >();
  test_100_segments_angles< CGAL::Exact_predicates_inexact_constructions_kernel >();
  test_100_segments_angles< CGAL::Exact_predicates_exact_constructions_kernel >();
  std::cout << "test_100_segments_angles: SUCCESS" << std::endl;
  return EXIT_SUCCESS;
}
