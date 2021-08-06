#include "include/utils.h"
#include "include/Saver.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Shape_regularization/regularize_segments.h>

namespace SR = CGAL::Shape_regularization;

template<class Traits>
void test_overloads() {

  using FT        = typename Traits::FT;
  using Point_2   = typename Traits::Point_2;
  using Segment_2 = typename Traits::Segment_2;
  using Segments  = std::vector<Segment_2>;

  using NQ = SR::Segments::Delaunay_neighbor_query_2<Traits, Segments>;
  using AR = SR::Segments::Angle_regularization_2<Traits, Segments>;
  using OR = SR::Segments::Offset_regularization_2<Traits, Segments>;
  using QP = CGAL::OSQP_quadratic_program_traits<FT>;

  Traits traits;
  Segments segments = {
    Segment_2(Point_2(0, 0), Point_2(1, 0)),
    Segment_2(Point_2(0, 1), Point_2(1, 1)),
  };
  assert(segments.size() == 2);

  // saver.export_segments(segments, "ol_input", 100);

  NQ neighbor_query(segments);
  AR angle_regularization(segments);
  OR offset_regularization(segments);
  QP quadratic_program;

  SR::Segments::regularize_segments(
    segments, neighbor_query, angle_regularization, quadratic_program, traits);
  SR::Segments::regularize_segments(
    segments, neighbor_query, angle_regularization, traits);
  SR::Segments::regularize_segments(
    segments, neighbor_query, angle_regularization);
  SR::Segments::regularize_segments(
    segments, traits);
  SR::Segments::regularize_segments(
    segments);

  SR::Segments::regularize_angles(
    segments, traits);
  SR::Segments::regularize_angles(
    segments);

  SR::Segments::regularize_offsets(
    segments, traits);
  SR::Segments::regularize_offsets(
    segments);
}

int main() {
  test_overloads< CGAL::Simple_cartesian<double> >();
  test_overloads< CGAL::Exact_predicates_inexact_constructions_kernel >();
  test_overloads< CGAL::Exact_predicates_exact_constructions_kernel >();
  std::cout << "test_overloads: SUCCESS" << std::endl;
  return EXIT_SUCCESS;
}
