#ifndef CGAL_ARRANGEMENT_ON_CURVE_1_TEST_UTILS_H
#define CGAL_ARRANGEMENT_ON_CURVE_1_TEST_UTILS_H

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <string>
#include <type_traits>
#include <vector>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arrangement_on_curve_1/Arrangement_on_curve_1.h>
#include <CGAL/Arrangement_on_curve_1/Line_traits_1.h>
#include <CGAL/Arrangement_on_curve_1/Segment_traits_1.h>
#include <CGAL/Arrangement_on_curve_1/Left_ray_traits_1.h>
#include <CGAL/Arrangement_on_curve_1/Right_ray_traits_1.h>
#include <CGAL/Arrangement_on_curve_1/Circle_traits_1.h>
#include <CGAL/Arrangement_on_curve_1/insert.h>
#include <CGAL/Arrangement_on_curve_1/locate.h>

namespace Aoc = CGAL::Arrangement_on_curve_1;

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using FT = Kernel::FT;
using Direction_2 = Kernel::Direction_2;

using Open_traits = Aoc::Line_traits_1<Kernel>;
using Closed_traits = Aoc::Segment_traits_1<Kernel>;
using Left_closed_traits = Aoc::Left_ray_traits_1<Kernel>;
using Right_closed_traits = Aoc::Right_ray_traits_1<Kernel>;
using Identified_traits = Aoc::Circle_traits_1<Kernel>;

// The arrangement on the geometry traits `Gt` with the default topology
// traits in list mode (`UseVector = false`) or vector mode (`UseVector = true`).
template <typename Gt, bool UseVector>
using Test_arrangement =
  Aoc::Arrangement_on_curve_1<Gt, typename Aoc::Default_topology_traits<Gt, void, void, UseVector>::type>;

// sorts points along the curve of `arr` using the comparison of its
// geometry traits; in an identified arrangement, the identification point is
// last (`Compare_x_1` must not be called on it).
template <typename Arr>
void sort_points(const Arr& arr, std::vector<typename Arr::Point_1>& pts) {
  const auto& gt = arr.geometry_traits_1();
  auto cmp = gt.compare_x_1_object();
  auto is_on_id = [&](const typename Arr::Point_1& p) {
    if constexpr (Arr::is_identified) return gt.is_on_x_identification_1_object()(p);
    else return false;
  };
  std::sort(pts.begin(), pts.end(), [&](const typename Arr::Point_1& a, const typename Arr::Point_1& b) {
    if (is_on_id(a)) return false;
    if (is_on_id(b)) return true;
    return cmp(a, b) == CGAL::SMALLER;
  });
}

// obtains the vertex of `arr` associated with `p`.
template <typename Arr>
typename Arr::Vertex_descriptor find_vertex(Arr& arr, const typename Arr::Point_1& p)
{ return std::get<0>(Aoc::locate(arr, p)); }

// calls `f` with every sequence of distinct indices in [0, n), that is, every
// permutation of every subset.
template <typename F>
void for_each_sequence(std::size_t n, F f) {
  for (unsigned mask = 0; mask < (1u << n); ++mask) {
    std::vector<int> seq;
    for (unsigned i = 0; i < n; ++i) if (mask & (1u << i)) seq.push_back(int(i));
    do f(seq);
    while (std::next_permutation(seq.begin(), seq.end()));
  }
}

// determines whether `arr` is a valid arrangement whose vertices are
// associated with the points in `expected` (in any order).  It checks the
// order of the vertices, the incidence relations, the number of edges, the
// boundary edges, the cached boundary and identification vertices, and
// that every point is located at its vertex.  If `arr` is not valid, `why`
// describes the first violation.
template <typename Arr>
bool is_valid(const Arr& arr, std::vector<typename Arr::Point_1> expected, std::string& why) {
  sort_points(arr, expected);
  const auto& tt = arr.topology_traits();
  auto pm = arr.vertex_point_map();
  auto lb = arr.left_boundary_edge();
  auto rb = arr.right_boundary_edge();
  if (arr.number_of_vertices() != expected.size()) { why = "number of vertices"; return false; }

  std::vector<typename Arr::Point_1> pts;
  std::size_t ne = 1;
  if constexpr (Arr::is_identified) {
    if (expected.empty()) {
      if ((arr.number_of_edges() != 1) || arr.has_left_vertex(lb) || arr.has_right_vertex(lb))
      { why = "empty arrangement"; return false; }
      if (tt.has_identification_vertex()) { why = "identification vertex of an empty arrangement"; return false; }
      if (rb != lb) { why = "right boundary edge of an empty arrangement"; return false; }
      return true;
    }
    auto e = lb;
    for (std::size_t i = 0; i < expected.size(); ++i) {
      if (! arr.has_left_vertex(e) || ! arr.has_right_vertex(e)) { why = "edge without an end vertex"; return false; }
      auto v = arr.right_vertex(e);
      if (arr.left_edge(v) != e) { why = "left edge"; return false; }
      pts.push_back(get(pm, v));
      e = arr.right_edge(v);
      if (arr.left_vertex(e) != v) { why = "left vertex"; return false; }
    }
    if (e != lb) { why = "the cycle is not closed"; return false; }
    ne = expected.size();
    const bool has_id = arr.geometry_traits_1().is_on_x_identification_1_object()(expected.back());
    if (tt.has_identification_vertex() != has_id) { why = "identification vertex flag"; return false; }
    if (has_id && (arr.left_vertex(lb) != tt.identification_vertex()))
    { why = "left boundary edge does not follow the identification vertex"; return false; }
    if (rb != (has_id ? arr.left_edge(tt.identification_vertex()) : lb))
    { why = "right boundary edge"; return false; }
  }
  else {
    const bool has_lbv = arr.has_left_vertex(lb);
    if constexpr (Arr::left_is_closed) {
      if (tt.has_left_boundary_vertex() != has_lbv) { why = "left boundary vertex flag"; return false; }
      if (has_lbv) {
        auto v = arr.left_vertex(lb);
        if ((tt.left_boundary_vertex() != v) || (arr.left_edge(v) != lb) || (arr.right_edge(v) != lb))
        { why = "left boundary vertex"; return false; }
        pts.push_back(get(pm, v));
      }
    }
    else if (has_lbv) { why = "left vertex of the left boundary edge on an open side"; return false; }

    auto e = lb;
    bool has_rbv = false;
    while (arr.has_right_vertex(e)) {
      auto v = arr.right_vertex(e);
      if (arr.left_edge(v) != e) { why = "left edge"; return false; }
      pts.push_back(get(pm, v));
      // Both edges of the right boundary vertex are its left edge.
      if (arr.right_edge(v) == e) {
        has_rbv = true;
        if constexpr (Arr::right_is_closed) {
          if (! tt.has_right_boundary_vertex() || (tt.right_boundary_vertex() != v))
          { why = "right boundary vertex"; return false; }
        }
        else { why = "right boundary vertex on an open side"; return false; }
        break;
      }
      e = arr.right_edge(v);
      if (! arr.has_left_vertex(e) || (arr.left_vertex(e) != v)) { why = "left vertex"; return false; }
      ++ne;
    }
    if (rb != e) { why = "right boundary edge"; return false; }
    if constexpr (Arr::right_is_closed) {
      if (! has_rbv && tt.has_right_boundary_vertex()) { why = "right boundary vertex flag"; return false; }
    }
  }

  if (arr.number_of_edges() != ne) { why = "number of edges"; return false; }
  if (pts != expected) { why = "order of the vertices"; return false; }
  for (const auto& p : expected) {
    auto loc = Aoc::locate(arr, p);
    if ((loc.index() != 0) || (get(pm, std::get<0>(loc)) != p)) { why = "locate"; return false; }
  }
  return true;
}

// reports the outcome of a test case, and counts it as a failure if needed.
inline bool report(bool ok, const std::string& name, const std::string& why, std::size_t& failures) {
  if (! ok) {
    ++failures;
    std::cerr << "FAIL " << name << ": " << why << std::endl;
  }
  return ok;
}

#endif
