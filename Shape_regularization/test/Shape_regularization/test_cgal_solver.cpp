#include "include/Saver.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Shape_regularization.h>

namespace SR = CGAL::Shape_regularization;

// This test is still here in case we will add a CGAL solver later, but for now
// it is removed and we do not aim adding it in the near future.

template<class Traits>
void test_cgal_solver() {

  using FT        = typename Traits::FT;
  using Point_2   = typename Traits::Point_2;
  using Segment_2 = typename Traits::Segment_2;
  using Segments  = std::vector<Segment_2>;
  using Indices   = std::vector<std::size_t>;
  using Saver     = SR::Tests::Saver<Traits>;

  using NQ = SR::Segments::Delaunay_neighbor_query_2<Traits, Segments>;
  using AR = SR::Segments::Angle_regularization_2<Traits, Segments>;
  using QP = CGAL::OSQP_quadratic_program_traits<FT>;

  using QP_AR = SR::QP_regularization<Traits, Segments, NQ, AR, QP>;

  Saver saver;
  Segments segments = {
    Segment_2(Point_2(0, 0), Point_2(1, 0)),
    Segment_2(Point_2(1, 0), Point_2(1, 1)),
    Segment_2(Point_2(1, 1), Point_2(0, 1)),
    Segment_2(Point_2(0, 1), Point_2(0, 0))
  };
  assert(segments.size() == 4);

  // saver.export_segments(segments, "cgal_input", 100);

  const FT max_angle_2 = FT(5);
  NQ neighbor_query(segments);
  AR angle_regularization(
    segments, CGAL::parameters::maximum_angle(max_angle_2));

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

  // saver.export_segments(segments, "cgal_angles", 100);

  assert(segments.size() == 4);
  assert(parallel_groups.size() == 2);
  assert(orthogonal_groups.size() == 1);
  assert(num_segments_angles == 4);
}

int main() {
  test_cgal_solver< CGAL::Simple_cartesian<double> >();
  test_cgal_solver< CGAL::Exact_predicates_inexact_constructions_kernel >();
  test_cgal_solver< CGAL::Exact_predicates_exact_constructions_kernel >();
  std::cout << "test_cgal_solver: SUCCESS" << std::endl;
  return EXIT_SUCCESS;
}
