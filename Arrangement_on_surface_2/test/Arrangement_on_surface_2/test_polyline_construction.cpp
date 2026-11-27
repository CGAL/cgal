// Test the construction of x-monotone polylines from points by the polyline
// traits, that is, the construction from a range of points and the appending
// of a point at the back or at the front of an existing x-monotone polyline.
// The subcurves must be constructed by the subcurve traits, as the subcurve
// type does not necessarily have a constructor from two points (e.g., a
// geodesic arc on the sphere). Every construction, of polylines directed left
// to right and right to left, is compared with the construction of the same
// polyline from a range of subcurves, which are constructed by the subcurve
// traits.

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_geodesic_arc_on_sphere_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Segment_traits = CGAL::Arr_segment_traits_2<Kernel>;
using Geodesic_traits = CGAL::Arr_geodesic_arc_on_sphere_traits_2<Kernel>;
using Segment_polyline_traits = CGAL::Arr_polyline_traits_2<Segment_traits>;
using Geodesic_polyline_traits = CGAL::Arr_polyline_traits_2<Geodesic_traits>;

/* Check whether two x-monotone polylines consist of the same subcurves in the
 * same order and in the same directions.
 */
template <typename Gt>
bool equal(const Gt& traits, const typename Gt::X_monotone_curve_2& xcv1,
           const typename Gt::X_monotone_curve_2& xcv2) {
  if (xcv1.number_of_subcurves() != xcv2.number_of_subcurves()) return false;
  const auto& sub_traits = *traits.subcurve_traits_2();
  auto equal_2 = sub_traits.equal_2_object();
  auto cmp_endpts = sub_traits.compare_endpoints_xy_2_object();
  for (std::size_t i = 0; i < xcv1.number_of_subcurves(); ++i) {
    if (! equal_2(xcv1[i], xcv2[i])) return false;
    if (cmp_endpts(xcv1[i]) != cmp_endpts(xcv2[i])) return false;
  }
  return true;
}

/* Construct the x-monotone polyline through the given points from the
 * subcurves constructed by the subcurve traits.
 */
template <typename Gt>
typename Gt::X_monotone_curve_2 expected_polyline(const Gt& traits, const std::vector<typename Gt::Point_2>& points) {
  auto ctr_sub = traits.subcurve_traits_2()->construct_x_monotone_curve_2_object();
  std::vector<typename Gt::X_monotone_subcurve_2> subcurves;
  for (std::size_t i = 0; i + 1 < points.size(); ++i) subcurves.push_back(ctr_sub(points[i], points[i + 1]));
  return traits.construct_x_monotone_curve_2_object()(subcurves.begin(), subcurves.end());
}

/* Test the constructions of the x-monotone polyline through the given points,
 * which must be sorted lexicographically (in either direction) and contain at
 * least three points.
 */
template <typename Gt>
bool test_constructions(const Gt& traits, const std::vector<typename Gt::Point_2>& points, const std::string& name) {
  using Point = typename Gt::Point_2;
  auto expected = expected_polyline(traits, points);
  auto ctr_xcv = traits.construct_x_monotone_curve_2_object();
  bool ok = true;

  // Construct from a range of points.
  auto xcv = ctr_xcv(points.begin(), points.end());
  if (! equal(traits, xcv, expected)) {
    std::cerr << "Error: " << name << ": construction from a range of points" << std::endl;
    ok = false;
  }

  // Append the last point at the back of the polyline through the other points.
  std::vector<Point> init(points.begin(), points.end() - 1);
  auto back = ctr_xcv(init.begin(), init.end());
  Point last = points.back();
  traits.push_back_2_object()(back, last);
  if (! equal(traits, back, expected)) {
    std::cerr << "Error: " << name << ": appending a point at the back" << std::endl;
    ok = false;
  }

  // Append the first point at the front of the polyline through the other points.
  std::vector<Point> tail(points.begin() + 1, points.end());
  auto front = ctr_xcv(tail.begin(), tail.end());
  Point first = points.front();
  traits.push_front_2_object()(front, first);
  if (! equal(traits, front, expected)) {
    std::cerr << "Error: " << name << ": appending a point at the front" << std::endl;
    ok = false;
  }

  return ok;
}

/* Test the constructions of x-monotone polylines through the given points,
 * which must be sorted lexicographically in increasing order, directed left to
 * right and right to left.
 */
template <typename Gt>
bool test_both_directions(const Gt& traits, std::vector<typename Gt::Point_2> points, const std::string& name) {
  bool ok = test_constructions(traits, points, name + " directed left to right");
  std::reverse(points.begin(), points.end());
  ok &= test_constructions(traits, points, name + " directed right to left");
  return ok;
}

int main() {
  bool ok = true;

  // Segments
  using Segment_point = Segment_traits::Point_2;
  Segment_polyline_traits segment_traits;
  std::vector<Segment_point> segment_points = {Segment_point(0, 0), Segment_point(1, 2), Segment_point(3, 1),
                                               Segment_point(4, 3)};
  ok &= test_both_directions(segment_traits, segment_points, "segment polyline");

  // Geodesic arcs on a small circle parallel to the XY-plane; the azimuth
  // increases.
  Geodesic_polyline_traits geodesic_traits;
  auto ctr_p = geodesic_traits.subcurve_traits_2()->construct_point_2_object();
  std::vector<Geodesic_traits::Point_2> geodesic_points;
  const double z = 0.5, r = std::sqrt(1.0 - z * z);
  for (double t : {0.2, 0.7, 1.2, 1.7}) geodesic_points.push_back(ctr_p(r * std::cos(t), r * std::sin(t), z));
  ok &= test_both_directions(geodesic_traits, geodesic_points, "geodesic polyline");

  if (! ok) return EXIT_FAILURE;
  std::cout << "Passed" << std::endl;
  return EXIT_SUCCESS;
}
