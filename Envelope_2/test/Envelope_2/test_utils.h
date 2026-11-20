#ifndef CGAL_ENVELOPE_2_TEST_UTILS_H
#define CGAL_ENVELOPE_2_TEST_UTILS_H

#include <algorithm>
#include <iostream>
#include <list>
#include <set>
#include <utility>
#include <variant>
#include <vector>

#include <CGAL/Exact_rational.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_curve_data_traits_2.h>
#include <CGAL/Envelope_2/Envelope_diagram_1.h>
#include <CGAL/Envelope_2/envelope_2.h>

// leda_rational, or Gmpq, or Quotient<MP_float>
using NT = CGAL::Exact_rational;
using Kernel = CGAL::Cartesian<NT>;
using Segment_traits_2 = CGAL::Arr_segment_traits_2<Kernel>;
using Traits_2 = CGAL::Arr_curve_data_traits_2<Segment_traits_2, int>;
using Point_2 = Traits_2::Point_2;
using Segment_2 = Segment_traits_2::Curve_2;
using Curve_2 = Traits_2::Curve_2;
using Diagram_1 = CGAL::Envelope_2::Envelope_diagram_1<Traits_2>;
using Curve_list = std::list<Curve_2>;

/*! checks if a curve with the same associated data as the input
 * curve is in the given range.
 * \param begin The beginning of the range.
 * \param end Past-the-end iterator.
 * \param c The curve, the data of which we are searching.
 * \return True if we found a curve with the same data.
 */
template <typename I>
bool find_curve(I begin, I end, const Curve_2& c) {
  while (begin != end) {
    if (begin->data() == c.data()) return true;
    ++begin;
  }
  return false;
}

/*! checks a bounded edge of the envelope of a given set of segments.
 * \param segs The segments.
 * \param diag The diagram.
 * \param e The edge.
 * \param is_lower Does the diagram represent the lower or the upper envelope.
 * \return Whether the edge is correct.
 */
template <typename Segments, typename Diagram>
bool check_edge(const Segments& segs, const Diagram& diag, typename Diagram::Edge_const_descriptor e, bool is_lower) {
  Kernel ker;
  auto midpoint = ker.construct_midpoint_2_object();
  auto min_ver = ker.construct_min_vertex_2_object();
  auto max_ver = ker.construct_max_vertex_2_object();
  auto comp_y_at_x = ker.compare_y_at_x_2_object();

  const Point_2& p_left = diag.point(diag.left_vertex(e));
  const Point_2& p_right = diag.point(diag.right_vertex(e));

  // Get the midpoint of the current edge.
  Point_2 p_mid = midpoint(p_left, p_right);

  const auto& e_curves = diag.edge_curves(e);

  // Check that all associated curves are equal at p_mid:
  for (auto it = e_curves.begin(); it != e_curves.end(); ++it) {
    if (comp_y_at_x(p_mid, e_curves.front(), *it) != CGAL::EQUAL) {
      std::cerr << "The edge (" << p_left << ") -> (" << p_right << ") is associated with "
                << "both of the following curves which are not equal: "
                << " [" << e_curves.front() << "], [" << *it << "]" << std::endl;
      return false;
    }
  }

  // Go over all segments.
  for (const auto& seg : segs) {
    // Check if p_mid lies in the x-range of the current segment.
    CGAL::Comparison_result res1 = CGAL::compare_x(p_mid, min_ver(seg));
    CGAL::Comparison_result res2 = CGAL::compare_x(p_mid, max_ver(seg));
    if (res1 == CGAL::SMALLER || res2 == CGAL::LARGER) continue;

    if (e_curves.empty()) {
      // If p_mid is in the x-range of the given segment, the current edge
      // cannot be empty.
      std::cerr << "The edge (" << p_left << ") -> (" << p_right << ") is empty, "
                << "but the segment [" << seg << "] is defined over it." << std::endl;
      return false;
    }

    // Compare the segment with the segment associated with the edge.
    CGAL::Comparison_result y_res = comp_y_at_x(p_mid, seg, e_curves.front());

    // Check for violations of the diagram properties.
    if ((is_lower && y_res == CGAL::SMALLER) || (! is_lower && y_res == CGAL::LARGER)) {
      std::cerr << "The edge (" << p_left << ") -> (" << p_right << ") is associated with the segment ["
                << e_curves.front() << "], but the segment [" << seg << "] violates the envelope properties."
                << std::endl;
      return false;
    }

    // If it is equal, the segment should be in the diagram.
    if (y_res == CGAL::EQUAL && find_curve(e_curves.begin(), e_curves.end(), seg) == false) {
      std::cerr << "The edge (" << p_left << ") -> (" << p_right << ") is associated with the segment ["
                << e_curves.front() << "], but the segment [" << seg << "] is not in the list of its curves."
                << std::endl;
      return false;
    }
  }
  return true;
}

/*! checks a vertex of the envelope of a given set of segments.
 * \param segs The segments.
 * \param diag The diagram.
 * \param v The vertex.
 * \param is_lower Does the diagram represent the lower or the upper envelope.
 * \return Whether the vertex is correct.
 */
template <typename Segments, typename Diagram>
bool check_vertex(const Segments& segs, const Diagram& diag, typename Diagram::Vertex_const_descriptor v,
                  bool is_lower) {
  Kernel ker;
  auto min_ver = ker.construct_min_vertex_2_object();
  auto max_ver = ker.construct_max_vertex_2_object();
  auto comp_y_at_x = ker.compare_y_at_x_2_object();

  const Point_2& p_v = diag.point(v);
  const auto& v_curves = diag.vertex_curves(v);

  // Check that all associated curves are equal at the vertex:
  for (auto it = v_curves.begin(); it != v_curves.end(); ++it) {
    if (comp_y_at_x(p_v, *it) != CGAL::EQUAL) {
      std::cerr << "The vertex (" << p_v << ") is associated with "
                << " [" << *it << "] but they are not equal at the vertex." << std::endl;
      return false;
    }
  }

  // Go over all segments.
  for (const auto& seg : segs) {
    // Check if v lies in the x-range of the current segment.
    CGAL::Comparison_result res1 = CGAL::compare_x(p_v, min_ver(seg));
    CGAL::Comparison_result res2 = CGAL::compare_x(p_v, max_ver(seg));
    if (res1 == CGAL::SMALLER || res2 == CGAL::LARGER) continue;

    if (v_curves.empty()) {
      // If v is in the x-range of the given segment, the current vertex
      // cannot be empty.
      std::cerr << "The vertex (" << p_v << ") is empty, " << "but the segment [" << seg << "] is defined over it."
                << std::endl;
      return false;
    }

    // Compare the segment with the segment associated with the vertex.
    CGAL::Comparison_result y_res = comp_y_at_x(p_v, seg);

    // Check for violations of the diagram properties.
    if ((is_lower && y_res == CGAL::LARGER) || (! is_lower && y_res == CGAL::SMALLER)) {
      std::cerr << "The vertex (" << p_v << ") and the segment [" << seg << "] violate the envelope properties."
                << std::endl;
      return false;
    }

    // If it is equal, the segment should be in the diagram.
    if (y_res == CGAL::EQUAL && find_curve(v_curves.begin(), v_curves.end(), seg) == false) {
      std::cerr << "The vertex (" << p_v << ") does not contain the segment ["
                << seg << "] in the list of its curves, but they have" << " equal values." << std::endl;
      return false;
    }
  }
  return true;
}

/*! checks the envelope of a given set of segments: every vertex, every
 * bounded edge, and the two unbounded edges.
 * \param segs The segments.
 * \param diag The diagram.
 * \param is_lower Does the diagram represent the lower or the upper envelope.
 * \return Whether the diagram structure is correct.
 */
template <typename Segments, typename Diagram>
bool check_envelope(const Segments& segs, const Diagram& diag, bool is_lower) {
  auto e = diag.leftmost();

  // If the diagram has no vertices, the segment set must be empty.
  if (e == diag.rightmost()) {
    if (! segs.empty()) {
      std::cerr << "The diagram is empty, but the segment set is not." << std::endl;
      return false;
    }
    return diag.empty_edge_curves(e);
  }

  // The segments are bounded, so the unbounded edges are not associated with curves.
  if (! diag.empty_edge_curves(diag.leftmost()) || ! diag.empty_edge_curves(diag.rightmost())) {
    std::cerr << "An unbounded edge is associated with curves." << std::endl;
    return false;
  }

  for (;;) {
    auto v = diag.right_vertex(e);
    if (! check_vertex(segs, diag, v, is_lower)) return false;
    e = diag.right_edge(v);
    if (e == diag.rightmost()) return true;
    if (! check_edge(segs, diag, e, is_lower)) return false;
  }
}

/*! obtains the envelope value of a segment at a given \f$x\f$-coordinate.
 * For a vertical segment, this is its lower (upper) endpoint.
 * \param seg The segment.
 * \param x The \f$x\f$-coordinate.
 * \param is_lower Is the envelope the lower or the upper envelope.
 * \param y Output: The value.
 * \return Whether `x` lies in the \f$x\f$-range of `seg`.
 */
inline bool value_at(const Curve_2& seg, const NT& x, bool is_lower, NT& y) {
  const Point_2& p = seg.source();
  const Point_2& q = seg.target();
  if ((x < (std::min)(p.x(), q.x())) || (x > (std::max)(p.x(), q.x()))) return false;
  if (p.x() == q.x()) y = is_lower ? (std::min)(p.y(), q.y()) : (std::max)(p.y(), q.y());
  else y = p.y() + (x - p.x()) * (q.y() - p.y()) / (q.x() - p.x());
  return true;
}

/*! computes the envelope of a given set of segments at a given
 * \f$x\f$-coordinate by brute force.
 * \param segs The segments.
 * \param x The \f$x\f$-coordinate.
 * \param is_lower Is the envelope the lower or the upper envelope.
 * \param with_vertical Whether vertical segments participate.
 * \param y Output: The envelope value.
 * \param ids Output: The data of the segments that attain it.
 * \return Whether any segment is defined at `x`.
 */
template <typename Segments>
bool envelope_at(const Segments& segs, const NT& x, bool is_lower, bool with_vertical, NT& y, std::set<int>& ids) {
  bool found = false;
  ids.clear();
  for (const auto& seg : segs) {
    if (! with_vertical && (seg.source().x() == seg.target().x())) continue;
    NT y_seg;
    if (! value_at(seg, x, is_lower, y_seg)) continue;
    if (! found || (is_lower ? (y_seg < y) : (y_seg > y))) {
      y = y_seg;
      ids.clear();
      found = true;
    }
    if (y_seg == y) ids.insert(seg.data());
  }
  return found;
}

/*! checks the envelope of a given set of segments against a brute-force
 * reference.  The candidate \f$x\f$-coordinates are those of all the
 * endpoints and of all the pairwise intersections.  At every candidate and at
 * the midpoint between every two consecutive candidates, the envelope value
 * and the set of segments that attain it are computed directly and compared
 * with the diagram feature at that \f$x\f$-coordinate.  The vertices must
 * lie at distinct candidates, in increasing \f$x\f$-order; a vertex must be at
 * the envelope and be associated with exactly the segments that attain it;
 * an edge must be associated with exactly the non-vertical segments that
 * attain the envelope, and a candidate in the interior of an edge must not
 * require a vertex.  The cost is quadratic in the number of segments.
 * \param segs The segments.
 * \param diag The diagram.
 * \param is_lower Does the diagram represent the lower or the upper envelope.
 * \return Whether the diagram matches the reference.
 */
template <typename Segments, typename Diagram>
bool check_envelope_reference(const Segments& segs, const Diagram& diag, bool is_lower) {
  using Kernel_segment_2 = Kernel::Segment_2;

  // Collect the candidate x-coordinates.
  std::set<NT> xs;
  std::vector<Kernel_segment_2> ksegs;
  for (const auto& seg : segs) {
    xs.insert(seg.source().x());
    xs.insert(seg.target().x());
    ksegs.push_back(Kernel_segment_2(seg.source(), seg.target()));
  }
  for (std::size_t i = 0; i < ksegs.size(); ++i) {
    for (std::size_t j = i + 1; j < ksegs.size(); ++j) {
      auto res = CGAL::intersection(ksegs[i], ksegs[j]);
      if (! res) continue;
      if (const auto* p = std::get_if<Point_2>(&*res)) xs.insert(p->x());
      else if (const auto* s = std::get_if<Kernel_segment_2>(&*res)) {
        xs.insert(s->source().x());
        xs.insert(s->target().x());
      }
    }
  }

  // The vertices lie at distinct candidates, in increasing x-order.
  {
    auto e = diag.leftmost();
    const NT* prev_x = nullptr;
    while (e != diag.rightmost()) {
      auto v = diag.right_vertex(e);
      const NT& x = diag.point(v).x();
      if ((xs.count(x) == 0) || (prev_x && (*prev_x >= x))) {
        std::cerr << "The vertex (" << diag.point(v) << ") is not at a distinct candidate x-coordinate." << std::endl;
        return false;
      }
      prev_x = &x;
      e = diag.right_edge(v);
    }
  }

  // The samples, in increasing x-order; the second element indicates a candidate.
  std::vector<NT> cands(xs.begin(), xs.end());
  std::vector<std::pair<NT, bool>> samples;
  for (std::size_t i = 0; i < cands.size(); ++i) {
    samples.push_back(std::make_pair(cands[i], true));
    if (i + 1 < cands.size()) samples.push_back(std::make_pair((cands[i] + cands[i + 1]) / 2, false));
  }

  // Walk along the diagram and the samples in parallel.
  auto e = diag.leftmost();
  for (const auto& [x, is_cand] : samples) {
    // Advance to the feature at x: the vertex at x, or the edge that contains x.
    bool on_vertex = false;
    typename Diagram::Vertex_const_descriptor v{};
    while (e != diag.rightmost()) {
      auto w = diag.right_vertex(e);
      CGAL::Comparison_result res = CGAL::compare(diag.point(w).x(), x);
      if (res == CGAL::LARGER) break;
      if (res == CGAL::EQUAL) {
        on_vertex = true;
        v = w;
        break;
      }
      e = diag.right_edge(w);
    }

    NT y;
    std::set<int> ids;
    const bool defined = envelope_at(segs, x, is_lower, on_vertex || is_cand, y, ids);
    if (on_vertex) {
      std::set<int> v_ids;
      for (const auto& cv : diag.vertex_curves(v)) v_ids.insert(cv.data());
      if (! defined || (diag.point(v).y() != y) || (v_ids != ids)) {
        std::cerr << "The vertex (" << diag.point(v) << ") does not match the reference envelope." << std::endl;
        return false;
      }
      continue;
    }

    std::set<int> e_ids;
    for (const auto& cv : diag.edge_curves(e)) e_ids.insert(cv.data());
    if (is_cand) {
      // A candidate in the interior of an edge: the non-vertical segments that
      // attain the envelope there must be those of the edge, and no other
      // segment (e.g., a vertical one or an endpoint) may attain it.
      NT y_e;
      std::set<int> ids_e;
      const bool defined_e = envelope_at(segs, x, is_lower, false, y_e, ids_e);
      if ((defined != defined_e) || (defined && ((ids != ids_e) || (y != y_e)))) {
        std::cerr << "A vertex is missing at x = " << x << "." << std::endl;
        return false;
      }
      if ((defined_e != ! e_ids.empty()) || (defined_e && (ids_e != e_ids))) {
        std::cerr << "The edge that contains x = " << x << " does not match the reference envelope." << std::endl;
        return false;
      }
      continue;
    }
    if ((defined != ! e_ids.empty()) || (defined && (e_ids != ids))) {
      std::cerr << "The edge that contains x = " << x << " does not match the reference envelope." << std::endl;
      return false;
    }
  }
  return true;
}

#endif
