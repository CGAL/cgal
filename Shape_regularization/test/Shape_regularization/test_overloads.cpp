#include "include/utils.h"
#include "include/Saver.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Shape_regularization/regularize_segments.h>
#include <CGAL/Shape_regularization/regularize_contours.h>

namespace SR = CGAL::Shape_regularization;

template<class Traits>
void test_overloads() {

  using FT        = typename Traits::FT;
  using Point_2   = typename Traits::Point_2;
  using Segment_2 = typename Traits::Segment_2;

  using Points   = std::vector<Point_2>;
  using Segments = std::vector<Segment_2>;
  using Indices  = std::vector<std::size_t>;

  using NQ = SR::Segments::Delaunay_neighbor_query_2<Traits, Segments>;
  using AR = SR::Segments::Angle_regularization_2<Traits, Segments>;
  using OR = SR::Segments::Offset_regularization_2<Traits, Segments>;
  using QP = CGAL::OSQP_quadratic_program_traits<FT>;
  using CD = SR::Contours::Longest_direction_2<Traits, Points>;

  Points points = {
    Point_2(0, 0), Point_2(1, 0),
    Point_2(0, 1), Point_2(1, 1),
  };
  assert(points.size() == 4);

  Segments segments = {
    Segment_2(Point_2(0, 0), Point_2(1, 0)),
    Segment_2(Point_2(0, 1), Point_2(1, 1)),
  };
  assert(segments.size() == 2);

  std::vector<Indices> indices;
  std::vector<Segment_2> unique;
  std::vector<Point_2> contour;

  // saver.export_segments(segments, "ol_input", 100);

  NQ neighbor_query(segments);
  AR angle_regularization(segments);
  OR offset_regularization(segments);
  QP quadratic_program;

  CD cdirections(points, true);
  CD odirections(points, false);

  SR::Segments::regularize_segments(
    segments, neighbor_query, angle_regularization, quadratic_program, CGAL::parameters::default_values());
  SR::Segments::regularize_segments(
    segments, neighbor_query, angle_regularization, quadratic_program);
  SR::Segments::regularize_segments(
    segments, neighbor_query, angle_regularization);
  SR::Segments::regularize_segments(
    segments);

  SR::Segments::regularize_angles(
    segments);
  SR::Segments::regularize_offsets(
    segments);

  SR::Segments::parallel_groups(
    segments, std::back_inserter(indices), CGAL::parameters::default_values());
  SR::Segments::parallel_groups(
    segments, std::back_inserter(indices));

  SR::Segments::collinear_groups(
    segments, std::back_inserter(indices), CGAL::parameters::default_values());
  SR::Segments::collinear_groups(
    segments, std::back_inserter(indices));

  SR::Segments::orthogonal_groups(
    segments, std::back_inserter(indices), CGAL::parameters::default_values());
  SR::Segments::orthogonal_groups(
    segments, std::back_inserter(indices));

  SR::Segments::unique_segments(
    segments, std::back_inserter(unique), CGAL::parameters::default_values());
  SR::Segments::unique_segments(
    segments, std::back_inserter(unique));

  SR::Contours::regularize_closed_contour(
    points, cdirections, std::back_inserter(contour), CGAL::parameters::default_values());
  SR::Contours::regularize_closed_contour(
    points, cdirections, std::back_inserter(contour));
  SR::Contours::regularize_closed_contour(
    points, std::back_inserter(contour));

  SR::Contours::regularize_open_contour(
    points, odirections, std::back_inserter(contour), CGAL::parameters::default_values());
  SR::Contours::regularize_open_contour(
    points, odirections, std::back_inserter(contour));
  SR::Contours::regularize_open_contour(
    points, std::back_inserter(contour));
}

int main() {
  test_overloads< CGAL::Simple_cartesian<double> >();
  test_overloads< CGAL::Exact_predicates_inexact_constructions_kernel >();
  test_overloads< CGAL::Exact_predicates_exact_constructions_kernel >();
  std::cout << "test_overloads: SUCCESS" << std::endl;
  return EXIT_SUCCESS;
}
