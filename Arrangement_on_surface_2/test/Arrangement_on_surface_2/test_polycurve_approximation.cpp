// Test the approximation of x-monotone polycurves by the polycurve traits
// (e.g., Arr_polyline_traits_2 and Arr_polycurve_traits_2). The approximation
// of a polycurve, in either direction, must be the concatenation of the
// approximations of its subcurves in the same direction, where the common
// endpoint of every two consecutive subcurves appears once. Also, test that
// the polycurve traits provides the approximation only to the extent that the
// subcurve traits does: not at all if the subcurve traits has no Approximate_2
// functor, and only of points if the Approximate_2 functor of the subcurve
// traits approximates only points.

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_circular_kernel_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_circular_arc_traits_2.h>
#include <CGAL/Arr_geodesic_arc_on_sphere_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arr_polycurve_traits_2.h>

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Segment_traits = CGAL::Arr_segment_traits_2<Kernel>;
using Geodesic_traits = CGAL::Arr_geodesic_arc_on_sphere_traits_2<Kernel>;
using Circular_arc_traits = CGAL::Arr_circular_arc_traits_2<CGAL::Exact_circular_kernel_2>;

using Segment_polyline_traits = CGAL::Arr_polyline_traits_2<Segment_traits>;
using Segment_polycurve_traits = CGAL::Arr_polycurve_traits_2<Segment_traits>;
using Geodesic_polyline_traits = CGAL::Arr_polyline_traits_2<Geodesic_traits>;
using Geodesic_polycurve_traits = CGAL::Arr_polycurve_traits_2<Geodesic_traits>;
using Nested_geodesic_polycurve_traits = CGAL::Arr_polycurve_traits_2<Geodesic_polycurve_traits>;

/* A segment traits whose Approximate_2 functor approximates only points. */
class Point_approximation_segment_traits : public Segment_traits {
public:
  class Approximate_2 {
  public:
    using Approximate_number_type = double;
    using Approximate_kernel = CGAL::Simple_cartesian<double>;
    using Approximate_point_2 = Approximate_kernel::Point_2;

    Approximate_number_type operator()(const Point_2& p, int i) const
    { return CGAL::to_double((i == 0) ? p.x() : p.y()); }

    Approximate_point_2 operator()(const Point_2& p) const
    { return Approximate_point_2(operator()(p, 0), operator()(p, 1)); }
  };

  Approximate_2 approximate_2_object() const { return Approximate_2(); }
};

using Point_approximation_polycurve_traits = CGAL::Arr_polycurve_traits_2<Point_approximation_segment_traits>;

// Detect whether the Approximate_2 functor of Gt approximates x-monotone curves.
template <typename Gt, typename = void>
struct has_approximate_xcv : std::false_type {};

template <typename Gt>
struct has_approximate_xcv
<Gt, std::void_t<decltype(std::declval<const typename Gt::Approximate_2&>()
                          (std::declval<const typename Gt::X_monotone_curve_2&>(), std::declval<double>(),
                           std::declval<std::back_insert_iterator<std::vector<typename Gt::Approximate_point_2>>>(),
                           std::declval<bool>()))>> :
  std::true_type
{};

// A polycurve of subcurves without Approximate_2 has no Approximate_2.
static_assert(std::is_void_v<CGAL::Arr_polycurve_traits_2<Circular_arc_traits>::Approximate_2>);

// A polycurve of subcurves whose Approximate_2 approximates only points
// approximates only points.
static_assert(! std::is_void_v<Point_approximation_polycurve_traits::Approximate_2>);
static_assert(! has_approximate_xcv<Point_approximation_polycurve_traits>::value);

// Polycurves of subcurves whose Approximate_2 approximates curves are approximated.
static_assert(has_approximate_xcv<Segment_polyline_traits>::value);
static_assert(has_approximate_xcv<Segment_polycurve_traits>::value);
static_assert(has_approximate_xcv<Geodesic_polyline_traits>::value);
static_assert(has_approximate_xcv<Geodesic_polycurve_traits>::value);
static_assert(has_approximate_xcv<Nested_geodesic_polycurve_traits>::value);

const double error = 0.01;

/* Approximate the given x-monotone curve. */
template <typename Gt>
std::vector<typename Gt::Approximate_2::Approximate_point_2>
approximate(const Gt& traits, const typename Gt::X_monotone_curve_2& xcv, bool l2r) {
  std::vector<typename Gt::Approximate_2::Approximate_point_2> points;
  traits.approximate_2_object()(xcv, error, std::back_inserter(points), l2r);
  return points;
}

/* Test that the approximation of an x-monotone polycurve, in both directions,
 * is the concatenation of the approximations of its subcurves in the same
 * direction, where the common endpoint of every two consecutive subcurves
 * appears once.
 */
template <typename Gt>
bool test_approximation(const Gt& traits, const typename Gt::X_monotone_curve_2& xcv, const std::string& name) {
  using Approximate_point = typename Gt::Approximate_point_2;
  const auto& sub_traits = *traits.subcurve_traits_2();
  bool directed_right = (traits.compare_endpoints_xy_2_object()(xcv) == CGAL::SMALLER);
  bool ok = true;
  for (bool l2r : {true, false}) {
    auto points = approximate(traits, xcv, l2r);

    // The subcurves are stored in the direction of the polycurve.
    std::vector<typename Gt::X_monotone_subcurve_2> subcurves(xcv.subcurves_begin(), xcv.subcurves_end());
    if (directed_right != l2r) std::reverse(subcurves.begin(), subcurves.end());
    std::vector<Approximate_point> expected;
    for (const auto& sub : subcurves) {
      auto sub_points = approximate(sub_traits, sub, l2r);
      auto begin = sub_points.begin();
      if (! expected.empty()) ++begin;
      expected.insert(expected.end(), begin, sub_points.end());
    }

    if (points == expected) continue;
    std::cerr << "Error: " << name << ": approximation (l2r = " << l2r << ") has " << points.size()
              << " points, expected " << expected.size() << std::endl;
    ok = false;
  }
  return ok;
}

/* Test the approximation of polylines and polycurves of segments. */
bool test_segments() {
  using Point = Segment_traits::Point_2;
  std::vector<Point> points = {Point(0, 0), Point(1, 2), Point(3, 1), Point(4, 3)};
  bool ok = true;

  // Polylines directed left to right and right to left
  Segment_polyline_traits polyline_traits;
  auto ctr_polyline = polyline_traits.construct_x_monotone_curve_2_object();
  auto polyline = ctr_polyline(points.begin(), points.end());
  ok &= test_approximation(polyline_traits, polyline, "segment polyline directed left to right");
  ok &= test_approximation(polyline_traits, ctr_polyline(points.rbegin(), points.rend()),
                           "segment polyline directed right to left");

  // The approximation of a polyline of segments consists of its points.
  std::vector<Segment_polyline_traits::Approximate_point_2> expected;
  auto approx = polyline_traits.approximate_2_object();
  for (const auto& p : points) expected.push_back(approx(p));
  if (approximate(polyline_traits, polyline, true) != expected) {
    std::cerr << "Error: the approximation of a segment polyline differs from its points" << std::endl;
    ok = false;
  }
  std::reverse(expected.begin(), expected.end());
  if (approximate(polyline_traits, polyline, false) != expected) {
    std::cerr << "Error: the reverse approximation of a segment polyline differs from its points" << std::endl;
    ok = false;
  }

  // Polycurves directed left to right and right to left
  Segment_polycurve_traits polycurve_traits;
  auto ctr_seg = polycurve_traits.subcurve_traits_2()->construct_x_monotone_curve_2_object();
  std::vector<Segment_traits::X_monotone_curve_2> segs, rsegs;
  for (std::size_t i = 0; i + 1 < points.size(); ++i) segs.push_back(ctr_seg(points[i], points[i + 1]));
  for (std::size_t i = points.size() - 1; i > 0; --i) rsegs.push_back(ctr_seg(points[i], points[i - 1]));
  auto ctr_polycurve = polycurve_traits.construct_x_monotone_curve_2_object();
  ok &= test_approximation(polycurve_traits, ctr_polycurve(segs.begin(), segs.end()),
                           "segment polycurve directed left to right");
  ok &= test_approximation(polycurve_traits, ctr_polycurve(rsegs.begin(), rsegs.end()),
                           "segment polycurve directed right to left");

  // A polycurve whose subcurve traits approximates only points
  Point_approximation_polycurve_traits point_approximation_traits;
  auto ap = point_approximation_traits.approximate_2_object()(points[1]);
  if ((ap.x() != 1) || (ap.y() != 2)) {
    std::cerr << "Error: wrong approximation of a point by a polycurve traits" << std::endl;
    ok = false;
  }
  return ok;
}

/* Test the approximation of polylines and (nested) polycurves of geodesic arcs
 * on a small circle parallel to the XY-plane.
 */
bool test_geodesic_arcs() {
  Geodesic_polyline_traits polyline_traits;
  auto ctr_p = polyline_traits.subcurve_traits_2()->construct_point_2_object();
  std::vector<Geodesic_traits::Point_2> points;
  const double z = 0.5, r = std::sqrt(1.0 - z * z);
  for (double t : {0.2, 0.7, 1.2, 1.7, 2.2}) points.push_back(ctr_p(r * std::cos(t), r * std::sin(t), z));
  bool ok = true;

  // Polylines directed left to right and right to left; the azimuth increases.
  auto ctr_arc = polyline_traits.subcurve_traits_2()->construct_x_monotone_curve_2_object();
  std::vector<Geodesic_traits::X_monotone_curve_2> arcs, rarcs;
  for (std::size_t i = 0; i + 1 < points.size(); ++i) arcs.push_back(ctr_arc(points[i], points[i + 1]));
  for (std::size_t i = points.size() - 1; i > 0; --i) rarcs.push_back(ctr_arc(points[i], points[i - 1]));
  auto ctr_polyline = polyline_traits.construct_x_monotone_curve_2_object();
  auto polyline = ctr_polyline(arcs.begin(), arcs.end());
  auto rpolyline = ctr_polyline(rarcs.begin(), rarcs.end());
  ok &= test_approximation(polyline_traits, polyline, "geodesic polyline directed left to right");
  ok &= test_approximation(polyline_traits, rpolyline, "geodesic polyline directed right to left");

  // The approximation of a geodesic arc consists of more than its endpoints.
  if (approximate(polyline_traits, polyline, true).size() <= points.size()) {
    std::cerr << "Error: the approximation of a geodesic polyline has too few points" << std::endl;
    ok = false;
  }

  // Polycurves of polycurves, each consisting of two arcs, directed left to
  // right and right to left
  Nested_geodesic_polycurve_traits nested_traits;
  const auto& polycurve_traits = *nested_traits.subcurve_traits_2();
  auto ctr_polycurve = polycurve_traits.construct_x_monotone_curve_2_object();
  auto ctr_nested = nested_traits.construct_x_monotone_curve_2_object();
  std::vector<Geodesic_polycurve_traits::X_monotone_curve_2> polycurves, rpolycurves;
  for (std::size_t i = 0; i + 1 < arcs.size(); i += 2) {
    std::vector<Geodesic_traits::X_monotone_curve_2> pair = {arcs[i], arcs[i + 1]};
    polycurves.push_back(ctr_polycurve(pair.begin(), pair.end()));
  }
  for (std::size_t i = 0; i + 1 < rarcs.size(); i += 2) {
    std::vector<Geodesic_traits::X_monotone_curve_2> pair = {rarcs[i], rarcs[i + 1]};
    rpolycurves.push_back(ctr_polycurve(pair.begin(), pair.end()));
  }
  auto nested = ctr_nested(polycurves.begin(), polycurves.end());
  auto rnested = ctr_nested(rpolycurves.begin(), rpolycurves.end());
  ok &= test_approximation(nested_traits, nested, "nested geodesic polycurve directed left to right");
  ok &= test_approximation(nested_traits, rnested, "nested geodesic polycurve directed right to left");

  // The approximation of a nested polycurve is identical to the approximation
  // of the polyline through the same points.
  for (bool l2r : {true, false}) {
    if ((approximate(nested_traits, nested, l2r) == approximate(polyline_traits, polyline, l2r)) &&
        (approximate(nested_traits, rnested, l2r) == approximate(polyline_traits, rpolyline, l2r)))
      continue;
    std::cerr << "Error: the approximations of a nested geodesic polycurve and of a geodesic polyline (l2r = "
              << l2r << ") differ" << std::endl;
    ok = false;
  }
  return ok;
}

int main() {
  bool ok = true;
  ok &= test_segments();
  ok &= test_geodesic_arcs();
  if (! ok) return EXIT_FAILURE;
  std::cout << "Passed" << std::endl;
  return EXIT_SUCCESS;
}
