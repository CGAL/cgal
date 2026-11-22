// Tests the configurations of the envelope algorithm that
// test_envelope_segments does not cover: the three pool strategies, vertical
// segments that induce the envelope at a vertex, degenerate configurations,
// the reuse of a diagram, the
// overloads that accept a traits object, constant input ranges, unbounded
// curves, and the deprecated interface.

#define CGAL_NO_DEPRECATION_WARNINGS

#include <algorithm>
#include <array>
#include <cstdlib>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include <CGAL/Arr_linear_traits_2.h>

#include "test_utils.h"

#ifndef CGAL_NO_DEPRECATED_CODE
#include <CGAL/Env_default_diagram_1.h>
#include <CGAL/envelope_2.h>
#endif

using Segment_vector = std::vector<Curve_2>;

// checks the envelope of the segments `segs` both structurally and against
// the brute-force reference.
template <typename Segments, typename Diagram>
bool is_valid_envelope(const Segments& segs, const Diagram& diag, bool is_lower) {
  return check_envelope(segs, diag, is_lower) && check_envelope_reference(segs, diag, is_lower);
}

// obtains a text representation of a diagram: the points of its vertices and
// the sorted data of the curves associated with its features, from left to
// right.
template <typename Diagram>
std::string dump(const Diagram& diag) {
  auto labels = [](const auto& curves) {
    std::vector<int> data;
    for (const auto& cv : curves) data.push_back(cv.data());
    std::sort(data.begin(), data.end());
    std::ostringstream os;
    os << '[';
    for (int d : data) os << ' ' << d;
    os << " ]";
    return os.str();
  };
  std::ostringstream os;
  auto e = diag.leftmost();
  for (;;) {
    os << 'E' << labels(diag.edge_curves(e));
    if (e == diag.rightmost()) break;
    auto v = diag.right_vertex(e);
    os << " V(" << diag.point(v) << ')' << labels(diag.vertex_curves(v)) << ' ';
    e = diag.right_edge(v);
  }
  return os.str();
}

// reports the outcome of a test case, and counts it as a failure if needed.
void report(bool ok, const std::string& name, std::size_t& failures) {
  std::cout << (ok ? "PASS " : "FAIL ") << name << std::endl;
  if (! ok) ++failures;
}

// generates `n` segments with integer endpoints in [0, range]^2, a quarter of
// which are vertical.  A small range yields degenerate configurations: shared
// endpoints, overlaps, and intersections at the x-coordinates of endpoints.
Segment_vector random_segments(std::mt19937& rng, int n, int range) {
  std::uniform_int_distribution<int> dist(0, range);
  Segment_vector segs;
  while (static_cast<int>(segs.size()) < n) {
    Point_2 p(dist(rng), dist(rng)), q(dist(rng), dist(rng));
    if (segs.size() % 4 == 3) q = Point_2(p.x(), q.y());
    if (p == q) continue;
    segs.push_back(Curve_2(Segment_2(p, q), static_cast<int>(segs.size())));
  }
  return segs;
}

// computes the envelope of `segs` with the given pool strategy, validates
// it, and returns its text representation.
template <bool Pooled, bool ValueBasedPool>
std::string envelope(Segment_vector& segs, bool is_lower, std::size_t& invalid) {
  using Dac = CGAL::Envelope_2::Envelope_divide_and_conquer_2<Traits_2, Diagram_1, Pooled, ValueBasedPool>;
  Diagram_1 diag;
  Dac env(is_lower ? Dac::LOWER : Dac::UPPER);
  env.insert_x_monotone_curves(segs.begin(), segs.end(), diag);
  if (! is_valid_envelope(segs, diag, is_lower)) ++invalid;
  return dump(diag);
}

// The value-based pool, the pointer-based pool, and no pool must all yield
// the same valid diagrams, also for degenerate inputs.
void test_pools(std::size_t& failures) {
  std::mt19937 rng(7);
  std::size_t invalid = 0, mismatches = 0;
  for (int t = 0; t < 200; ++t) {
    auto segs = random_segments(rng, 1 + t % 40, 3 + t % 18);
    for (bool is_lower : {true, false}) {
      auto value_based = envelope<true, true>(segs, is_lower, invalid);
      auto pointer_based = envelope<true, false>(segs, is_lower, invalid);
      auto no_pool = envelope<false, false>(segs, is_lower, invalid);
      if ((value_based != pointer_based) || (value_based != no_pool)) ++mismatches;
    }
  }
  report((invalid == 0) && (mismatches == 0), "pool strategies", failures);
}

// A vertical segment that lies below (above) a vertex of the lower (upper)
// envelope induces the envelope at the x-coordinate of the vertex, so the
// vertex moves to the endpoint of the vertical segment.
void test_vertical(std::size_t& failures) {
  Segment_vector lower_segs{Curve_2(Segment_2(Point_2(17, 3), Point_2(20, 3)), 0),
                            Curve_2(Segment_2(Point_2(20, 10), Point_2(20, 2)), 1)};
  Segment_vector upper_segs{Curve_2(Segment_2(Point_2(8, 3), Point_2(20, 9)), 0),
                            Curve_2(Segment_2(Point_2(20, 16), Point_2(20, 15)), 1)};
  Diagram_1 min_diag, max_diag;
  CGAL::Envelope_2::lower_envelope_x_monotone_2(lower_segs.begin(), lower_segs.end(), min_diag);
  CGAL::Envelope_2::upper_envelope_x_monotone_2(upper_segs.begin(), upper_segs.end(), max_diag);
  bool ok = is_valid_envelope(lower_segs, min_diag, true) && is_valid_envelope(upper_segs, max_diag, false);
  report(ok, "vertical segments at vertices", failures);
}

// constructs segments from the coordinates of their endpoints.
Segment_vector segments(const std::vector<std::array<int, 4>>& coords) {
  Segment_vector segs;
  for (const auto& c : coords)
    segs.push_back(Curve_2(Segment_2(Point_2(c[0], c[1]), Point_2(c[2], c[3])), static_cast<int>(segs.size())));
  return segs;
}

// Degenerate configurations that used to yield invalid diagrams or to fail
// an assertion.
void test_degenerate(std::size_t& failures) {
  struct Case {
    const char* name;
    bool is_lower;
    Segment_vector segs;
  };
  const std::vector<Case> cases{
    // An intersection point above the reached vertex yielded a second vertex at its x-coordinate.
    {"intersection above the reached vertex", true, segments({{5, 3, 15, 13}, {4, 8, 9, 1}, {3, 13, 11, 5}})},
    // The same configuration failed an assertion.
    {"intersection at the x-coordinate of an endpoint", true, segments({{15, 2, 18, 5}, {20, 2, 17, 4}, {4, 11, 17, 1}})},
    // An intersection point below the next vertex of the upper envelope yielded a second vertex.
    {"intersection below the next vertex", false, segments({{2, 2, 0, 2}, {0, 3, 1, 2}, {3, 0, 1, 3}})},
    // Curves that overlap up to the x-coordinate of a vertex were associated with it.
    {"overlap up to a vertex", true, segments({{2, 1, 4, 1}, {4, 2, 3, 0}, {0, 1, 3, 1}})},
    // A vertical segment on overlapping curves was associated with only one of them.
    {"vertical segment on an overlap (lower)", true, segments({{2, 1, 5, 1}, {4, 1, 0, 1}, {3, 2, 3, 1}})},
    {"vertical segment on an overlap (upper)", false, segments({{2, 3, 0, 3}, {2, 3, 0, 3}, {1, 3, 1, 2}})}};
  for (const auto& c : cases) {
    Diagram_1 diag;
    if (c.is_lower) CGAL::Envelope_2::lower_envelope_x_monotone_2(c.segs.begin(), c.segs.end(), diag);
    else CGAL::Envelope_2::upper_envelope_x_monotone_2(c.segs.begin(), c.segs.end(), diag);
    report(is_valid_envelope(c.segs, diag, c.is_lower), c.name, failures);
  }
}

// A diagram that is reused must not retain its previous content, also when
// the new input consists of vertical segments only, or is empty.
void test_reuse(std::size_t& failures) {
  std::mt19937 rng(11);
  auto segs1 = random_segments(rng, 30, 1000);
  Segment_vector segs2{Curve_2(Segment_2(Point_2(5, 1), Point_2(5, 8)), 0),
                       Curve_2(Segment_2(Point_2(9, 2), Point_2(9, 3)), 1)};
  auto segs3 = random_segments(rng, 10, 1000);
  Segment_vector none;

  Diagram_1 diag;
  CGAL::Envelope_2::lower_envelope_x_monotone_2(segs1.begin(), segs1.end(), diag);
  bool ok = is_valid_envelope(segs1, diag, true);
  CGAL::Envelope_2::lower_envelope_x_monotone_2(segs2.begin(), segs2.end(), diag);
  ok = ok && is_valid_envelope(segs2, diag, true);
  CGAL::Envelope_2::upper_envelope_x_monotone_2(segs3.begin(), segs3.end(), diag);
  ok = ok && is_valid_envelope(segs3, diag, false);
  CGAL::Envelope_2::lower_envelope_x_monotone_2(none.begin(), none.end(), diag);
  ok = ok && (diag.leftmost() == diag.rightmost()) && diag.empty_edge_curves(diag.leftmost());
  report(ok, "diagram reuse", failures);
}

// The overloads that accept a traits object, and the diagrams that share it.
// Without a traits object, the traits of the diagram are used.
void test_traits(std::size_t& failures) {
  std::mt19937 rng(13);
  auto segs = random_segments(rng, 30, 1000);
  auto traits = std::make_shared<const Traits_2>();
  Diagram_1 min_diag(traits), max_diag(traits);
  CGAL::Envelope_2::lower_envelope_x_monotone_2(segs.begin(), segs.end(), min_diag, traits);
  CGAL::Envelope_2::upper_envelope_x_monotone_2(segs.begin(), segs.end(), max_diag, traits);
  Diagram_1 diag(traits);
  CGAL::Envelope_2::lower_envelope_x_monotone_2(segs.begin(), segs.end(), diag);
  bool ok = is_valid_envelope(segs, min_diag, true) && is_valid_envelope(segs, max_diag, false) &&
    (dump(diag) == dump(min_diag)) && (min_diag.shared_geometry_traits_2() == traits) &&
    (&min_diag.geometry_traits_2() == traits.get());
  report(ok, "traits overloads", failures);
}

// The input range may be constant.
void test_const_range(std::size_t& failures) {
  std::mt19937 rng(19);
  const Segment_vector segs = random_segments(rng, 30, 1000);
  const Curve_list curves(segs.begin(), segs.end());
  Diagram_1 min_diag, max_diag;
  CGAL::Envelope_2::lower_envelope_x_monotone_2(segs.begin(), segs.end(), min_diag);
  CGAL::Envelope_2::upper_envelope_2(curves.begin(), curves.end(), max_diag);
  bool ok = is_valid_envelope(segs, min_diag, true) && is_valid_envelope(segs, max_diag, false);
  report(ok, "constant input ranges", failures);
}

// The envelopes of unbounded curves: the lines y = x (0) and y = -x (1), the
// ray y = -5 for x >= 0 (2), and the segment y = -9 for -3 <= x <= -1 (3).
void test_unbounded(std::size_t& failures) {
  using Linear_traits_2 = CGAL::Arr_curve_data_traits_2<CGAL::Arr_linear_traits_2<Kernel>, int>;
  using Linear_curve_2 = Linear_traits_2::X_monotone_curve_2;
  using Linear_diagram_1 = CGAL::Envelope_2::Envelope_diagram_1<Linear_traits_2>;
  using Line_2 = Kernel::Line_2;
  using Ray_2 = Kernel::Ray_2;
  using Kernel_segment_2 = Kernel::Segment_2;

  std::vector<Linear_curve_2> curves{Linear_curve_2(Line_2(Point_2(0, 0), Point_2(1, 1)), 0),
                                     Linear_curve_2(Line_2(Point_2(0, 0), Point_2(1, -1)), 1),
                                     Linear_curve_2(Ray_2(Point_2(0, -5), Point_2(1, -5)), 2),
                                     Linear_curve_2(Kernel_segment_2(Point_2(-3, -9), Point_2(-1, -9)), 3)};
  Linear_diagram_1 min_diag, max_diag;
  CGAL::Envelope_2::lower_envelope_x_monotone_2(curves.begin(), curves.end(), min_diag);
  CGAL::Envelope_2::upper_envelope_x_monotone_2(curves.begin(), curves.end(), max_diag);

  const std::string min_expected =
    "E[ 0 ] V(-3 -9)[ 3 ] E[ 3 ] V(-1 -9)[ 3 ] E[ 0 ] V(0 -5)[ 2 ] E[ 2 ] V(5 -5)[ 1 2 ] E[ 1 ]";
  const std::string max_expected = "E[ 1 ] V(0 0)[ 0 1 ] E[ 0 ]";
  bool ok = true;
  if (dump(min_diag) != min_expected) {
    std::cerr << "lower envelope: " << dump(min_diag) << std::endl;
    ok = false;
  }
  if (dump(max_diag) != max_expected) {
    std::cerr << "upper envelope: " << dump(max_diag) << std::endl;
    ok = false;
  }
  report(ok, "unbounded curves", failures);
}

#ifndef CGAL_NO_DEPRECATED_CODE
// The deprecated diagram class and the deprecated functions in the CGAL namespace.
void test_deprecated(std::size_t& failures) {
  std::mt19937 rng(17);
  auto segs = random_segments(rng, 25, 1000);
  Curve_list curves(segs.begin(), segs.end());
  Traits_2 traits;
  CGAL::Env_default_diagram_1<Traits_2> min_diag, max_diag;
  Diagram_1 diag;
  CGAL::lower_envelope_x_monotone_2(segs.begin(), segs.end(), min_diag, traits);
  CGAL::upper_envelope_x_monotone_2(segs.begin(), segs.end(), max_diag);
  CGAL::lower_envelope_2(curves.begin(), curves.end(), diag);
  bool ok = is_valid_envelope(segs, min_diag, true) && is_valid_envelope(segs, max_diag, false) &&
    (dump(diag) == dump(min_diag));
  report(ok, "deprecated interface", failures);
}
#endif

int main() {
  std::size_t failures = 0;
  test_pools(failures);
  test_vertical(failures);
  test_degenerate(failures);
  test_reuse(failures);
  test_traits(failures);
  test_const_range(failures);
  test_unbounded(failures);
#ifndef CGAL_NO_DEPRECATED_CODE
  test_deprecated(failures);
#endif
  return (failures == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
