#include "include/utils.h"
#include "include/Saver.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Shape_regularization.h>

// Typedefs.
using Kernel      = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT          = typename Kernel::FT;
using Segment_2   = typename Kernel::Segment_2;
using Segments    = std::vector<Segment_2>;
using Indices     = std::vector<std::size_t>;
using Segment_map = CGAL::Identity_property_map<Segment_2>;

using Neighbor_query =
  CGAL::Shape_regularization::Segments::Delaunay_neighbor_query_2<Kernel, Segments, Segment_map>;
using Angle_regularization =
  CGAL::Shape_regularization::Segments::Angle_regularization_2<Kernel, Segments, Segment_map>;
using Offset_regularization =
  CGAL::Shape_regularization::Segments::Offset_regularization_2<Kernel, Segments, Segment_map>;
using Quadratic_program =
  CGAL::OSQP_quadratic_program_traits<FT>;

using Quadratic_angle_regularizer =
  CGAL::Shape_regularization::QP_regularization<
    Kernel, Segments, Neighbor_query, Angle_regularization, Quadratic_program>;
using Quadratic_offset_regularizer =
  CGAL::Shape_regularization::QP_regularization<
    Kernel, Segments, Neighbor_query, Offset_regularization, Quadratic_program>;

int main(int argc, char *argv[]) {

  // If we want to save the result in a file, we save it in a path.
  std::string path = "";
  if (argc > 1) path = argv[1];
  Saver<Kernel> saver;

  // Initialize 15 segments.
  std::vector<Segment_2> segments;
  create_example_15(segments);

  // We create three groups of segments:
  // outer, top and bottom rhombuses.
  std::vector<Indices> groups(3);
  groups[0] = {0, 1, 2, 3, 4, 5, 6}; // outer
  groups[1] = {7, 8, 9, 10};         // top rhombus
  groups[2] = {11, 12, 13, 14};      // bottom rhombus

  // Save input segments.
  if (path != "") {
    const std::string full_path = path + "regularize_15_segments_before";
    saver.export_eps_segments(segments, full_path, FT(100));
  }

  // Angle regularization.
  const FT max_angle_2 = FT(10);

  // Create qp solver, neighbor query, and angle-based regularization model.
  Quadratic_program qp_angles;
  Neighbor_query neighbor_query(segments);
  Angle_regularization angle_regularization(
    segments, CGAL::parameters::maximum_angle(max_angle_2));

  // Add each group of input segments.
  for (const auto& group : groups) {
    neighbor_query.add_group(group);
    angle_regularization.add_group(group);
  }

  // Regularize.
  Quadratic_angle_regularizer qp_angle_regularizer(
    segments, neighbor_query, angle_regularization, qp_angles);
  qp_angle_regularizer.regularize();

  std::cout << "* number of modified segments (angles) = " <<
    angle_regularization.number_of_modified_segments() << std::endl;

  // Offset regularization.
  const FT max_offset_2 = FT(1) / FT(10);

  // Get groups of parallel segments after angle regularization.
  std::vector<Indices> pgroups;
  angle_regularization.parallel_groups(
    std::back_inserter(pgroups));

  // Create qp solver and offset-based regularization model.
  Quadratic_program qp_offsets;
  Offset_regularization offset_regularization(
    segments, CGAL::parameters::maximum_offset(max_offset_2));

  // Add each group of parallel segments with at least 2 segments.
  neighbor_query.clear();
  for (const auto& pgroup : pgroups) {
    neighbor_query.add_group(pgroup);
    offset_regularization.add_group(pgroup);
  }

  // Regularize.
  Quadratic_offset_regularizer qp_offset_regularizer(
    segments, neighbor_query, offset_regularization, qp_offsets);
  qp_offset_regularizer.regularize();

  std::cout << "* number of modified segments (offsets) = " <<
    offset_regularization.number_of_modified_segments() << std::endl;

  // Save regularized segments.
  if (path != "") {
    const std::string full_path = path + "regularize_15_segments_after";
    saver.export_eps_segments(segments, full_path, FT(100));
  }
}
