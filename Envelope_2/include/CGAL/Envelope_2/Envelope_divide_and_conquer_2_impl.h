// Copyright (c) 2006  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Ron Wein   <wein@post.tau.ac.il>
//             Efi Fogel  <efifogel@gmail.com>

#ifndef CGAL_ENVELOPE_DIVIDE_AND_CONQUER_2_IMPL_H
#define CGAL_ENVELOPE_DIVIDE_AND_CONQUER_2_IMPL_H

#include <CGAL/license/Envelope_2.h>

/*! \file
 * Definitions of the functions of the Envelope_divide_and_conquer_2 class.
 *
 * All accesses to vertex and edge data (curve lists, connectivity) are made
 * through the diagram's public interface rather than through raw node-pointer
 * arrow operators, because the new Envelope_diagram_1 derives from
 * Arrangement_on_curve_1 and stores curve data in property maps.
 *
 * Navigation shortcuts (v_left, v_right, e_left, e_right, …) are static
 * helpers defined in Env_divide_and_conquer_2.h.
 */


#include <optional>
#include <algorithm>

namespace CGAL {

// ---------------------------------------------------------------------------
// Construct the lower/upper envelope of non-vertical curves.
//
template <typename Traits, typename Diagram>
void Envelope_divide_and_conquer_2<Traits,Diagram>::
_construct_envelope_non_vertical(Curve_pointer_iterator begin, Curve_pointer_iterator end, Envelope_diagram_1& out_d) {
  out_d.clear();

  if (begin == end) return;

  // Check if the range contains just a single curve.
  Curve_pointer_iterator iter = begin;
  ++iter;

  if (iter == end) {
    // Construct a singleton diagram, which matches a single curve.
    _construct_singleton_diagram(*(*begin), out_d);
  }
  else {
    // Divide the given range of curves into two.
    std::size_t size = std::distance(begin, end);
    Curve_pointer_iterator div_it = begin;
    std::advance(div_it, size / 2);

    // Construct the diagrams (envelopes) for the two sub-ranges recursively
    // and then merge the two diagrams to obtain the result.
    Envelope_diagram_1 d1(out_d.shared_geometry_traits_1());
    Envelope_diagram_1 d2(out_d.shared_geometry_traits_1());

    _construct_envelope_non_vertical(begin, div_it, d1);
    _construct_envelope_non_vertical(div_it, end, d2);

    _merge_envelopes(d1, d2, out_d);
  }
}

// ---------------------------------------------------------------------------
// Construct a singleton diagram, which matches a single curve.
//
template <typename Traits, typename Diagram>
void Envelope_divide_and_conquer_2<Traits,Diagram>::
_construct_singleton_diagram(const X_monotone_curve_2& cv, Envelope_diagram_1& out_d) {
  CGAL_assertion(out_d.leftmost() == out_d.rightmost());
  CGAL_assertion(out_d.empty_edge_curves(out_d.leftmost()));

  auto& topo = out_d.topology_traits();

  // Check if the given curve is bounded from the left and from the right.
  if (m_traits->parameter_space_in_x_2_object()(cv, ARR_MIN_END) != ARR_INTERIOR) {
    if (m_traits->parameter_space_in_x_2_object()(cv, ARR_MAX_END) != ARR_INTERIOR) {
      // The curve is defined over (-oo, oo), so its diagram contains
      // only a single edge.
      out_d.add_edge_curve(out_d.leftmost(), cv);
      return;
    }

    // The curve is defined over (-oo, x], where x is finite.
    // Create a vertex and associate it with the right endpoint of cv.
    CGAL_precondition(m_traits->parameter_space_in_y_2_object()(cv, ARR_MAX_END) == ARR_INTERIOR);

    Vertex_handle v = out_d.new_vertex(m_traits->construct_max_vertex_2_object()(cv));
    Edge_handle e_right = out_d.new_edge();

    out_d.add_vertex_curve(v, cv);
    topo.set_left_edge(v, out_d.leftmost());
    topo.set_right_edge(v, e_right);

    out_d.add_edge_curve(out_d.leftmost(), cv);
    topo.set_right_vertex(out_d.leftmost(), v);

    topo.set_left_vertex(e_right, v);
    out_d.set_rightmost(e_right);
    return;
  }

  if (m_traits->parameter_space_in_x_2_object()(cv, ARR_MAX_END) != ARR_INTERIOR) {
    // The curve is defined over [x, +oo), where x is finite.
    // Create a vertex and associate it with the left endpoint of cv.
    CGAL_precondition(m_traits->parameter_space_in_y_2_object()(cv, ARR_MIN_END) == ARR_INTERIOR);

    Vertex_handle v = out_d.new_vertex(m_traits->construct_min_vertex_2_object()(cv));
    Edge_handle e_left = out_d.new_edge();

    out_d.add_vertex_curve(v, cv);
    topo.set_left_edge(v, e_left);
    topo.set_right_edge(v, out_d.rightmost());

    out_d.add_edge_curve(out_d.rightmost(), cv);
    topo.set_left_vertex(out_d.rightmost(), v);

    topo.set_right_vertex(e_left, v);
    out_d.set_leftmost(e_left);
    return;
  }

  // If we reached here, the curve is defined over a bounded x-range: [x1, x2]
  CGAL_precondition(m_traits->parameter_space_in_y_2_object()(cv, ARR_MIN_END) == ARR_INTERIOR);
  CGAL_precondition(m_traits->parameter_space_in_y_2_object()(cv, ARR_MAX_END) == ARR_INTERIOR);

  Vertex_handle v1 = out_d.new_vertex(m_traits->construct_min_vertex_2_object()(cv));
  Vertex_handle v2 = out_d.new_vertex(m_traits->construct_max_vertex_2_object()(cv));

  Edge_handle e_left = out_d.new_edge();
  Edge_handle e_right = out_d.new_edge();
  Edge_handle e = out_d.leftmost();

  out_d.add_vertex_curve(v1, cv);
  topo.set_left_edge(v1, e_left);
  topo.set_right_edge(v1, e);

  out_d.add_vertex_curve(v2, cv);
  topo.set_left_edge(v2, e);
  topo.set_right_edge(v2, e_right);

  out_d.add_edge_curve(e, cv);
  topo.set_left_vertex(e, v1);
  topo.set_right_vertex(e, v2);

  topo.set_right_vertex(e_left, v1);
  topo.set_left_vertex(e_right, v2);

  out_d.set_leftmost(e_left);
  out_d.set_rightmost(e_right);
}

// ---------------------------------------------------------------------------
// Merge two minimization (or maximization) diagrams.
//
template <typename Traits, typename Diagram>
void Envelope_divide_and_conquer_2<Traits,Diagram>::
_merge_envelopes(const Envelope_diagram_1& d1, const Envelope_diagram_1& d2, Envelope_diagram_1& out_d) {
  Edge_const_handle e1 = d1.leftmost();
  bool is_leftmost1 = true;
  Vertex_const_handle v1{};
  Edge_const_handle e2 = d2.leftmost();
  bool is_leftmost2 = true;
  Vertex_const_handle v2{};
  Vertex_const_handle next_v{};
  bool next_exists = true;
  Comparison_result res_v = EQUAL;
  bool same_x = false;

  do {
    same_x = false;

    if (e1 == d1.rightmost()) {
      if (e2 == d2.rightmost()) next_exists = false;
      else {
        v2 = d2.right_vertex(e2);
        next_v = v2;
        res_v = LARGER;
      }
    }
    else if (e2 == d2.rightmost()) {
      v1 = d1.right_vertex(e1);
      next_v = v1;
      res_v = SMALLER;
    }
    else {
      v1 = d1.right_vertex(e1);
      v2 = d2.right_vertex(e2);
      res_v = _compare_vertices(d1, d2, v1, v2, same_x);
      next_v = (res_v == SMALLER) ? v1 : v2;
    }

    if (! d1.empty_edge_curves(e1) && ! d2.empty_edge_curves(e2))
      _merge_two_intervals(e1, is_leftmost1, e2, is_leftmost2, next_v, next_exists, res_v, d1, d2, out_d);
    else if (! d1.empty_edge_curves(e1) && d2.empty_edge_curves(e2))
      _merge_single_interval(e1, e2, next_v, next_exists, res_v, d1, d2, out_d);
    else if (d1.empty_edge_curves(e1) && ! d2.empty_edge_curves(e2))
      _merge_single_interval(e2, e1, next_v, next_exists, CGAL::opposite(res_v), d2, d1, out_d);
    else {
      // Both empty.
      if (next_exists) {
        const Point_2& p_next = (res_v == SMALLER) ? d1.point(v1) : d2.point(v2);
        Vertex_handle new_v = _append_vertex(out_d, p_next, e1, d1);
        switch(res_v) {
         case SMALLER: out_d.add_vertex_curves(new_v, d1.vertex_curves(v1).begin(), d1.vertex_curves(v1).end()); break;
         case LARGER: out_d.add_vertex_curves(new_v, d2.vertex_curves(v2).begin(), d2.vertex_curves(v2).end()); break;
         case EQUAL:
          out_d.add_vertex_curves(new_v, d1.vertex_curves(v1).begin(), d1.vertex_curves(v1).end());
          out_d.add_vertex_curves(new_v, d2.vertex_curves(v2).begin(), d2.vertex_curves(v2).end());
          break;
        }
      }
    }

    if (next_exists) {
      if (res_v == SMALLER) {
        e1 = d1.right_edge(v1);
        is_leftmost1 = false;

        if (same_x) {
          e2 = d2.right_edge(v2);
          is_leftmost2 = false;
        }
      }
      else if (res_v == LARGER) {
        e2 = d2.right_edge(v2);
        is_leftmost2 = false;

        if (same_x) {
          e1 = d1.right_edge(v1);
          is_leftmost1 = false;
        }
      }
      else {
        e1 = d1.right_edge(v1);
        is_leftmost1 = false;

        e2 = d2.right_edge(v2);
        is_leftmost2 = false;
      }
    }

  } while (next_exists);
}

// ---------------------------------------------------------------------------
// Compare two diagram vertices.
//
template <typename Traits, typename Diagram>
Comparison_result Envelope_divide_and_conquer_2<Traits,Diagram>::
_compare_vertices(const Envelope_diagram_1& d1, const Envelope_diagram_1& d2,
                  Vertex_const_handle v1, Vertex_const_handle v2, bool& same_x) const {
  Comparison_result res = m_traits->compare_x_2_object()(d1.point(v1), d2.point(v2));

  if (res != EQUAL) {
    same_x = false;
    return res;
  }

  same_x = true;
  res = m_traits->compare_xy_2_object()(d1.point(v1), d2.point(v2));
  return (m_env_type == UPPER) ? CGAL::opposite(res) : res;
}

// ---------------------------------------------------------------------------
// Deal with an interval which is non-empty in one diagram and empty in the other.
//
template <typename Traits, typename Diagram>
void Envelope_divide_and_conquer_2<Traits,Diagram>::
_merge_single_interval(Edge_const_handle e, Edge_const_handle other_edge, Vertex_const_handle v, bool v_exists,
                       Comparison_result origin_of_v, const Envelope_diagram_1& in_d, const Envelope_diagram_1& other_d,
                       Envelope_diagram_1& out_d) {
  if (! v_exists) {
    out_d.add_edge_curves(out_d.rightmost(), in_d.edge_curves(e).begin(), in_d.edge_curves(e).end());
    return;
  }

  if (origin_of_v == SMALLER) {
    Vertex_handle new_v = _append_vertex(out_d, in_d.point(v), e, in_d);
    out_d.add_vertex_curves(new_v, in_d.vertex_curves(v).begin(), in_d.vertex_curves(v).end());
    return;
  }

  if (origin_of_v == EQUAL) {
    Vertex_handle new_v = _append_vertex(out_d, in_d.point(v), e, in_d);
    Vertex_const_handle v1 = in_d.right_vertex(e);
    Vertex_const_handle v2 = other_d.right_vertex(other_edge);
    out_d.add_vertex_curves(new_v, in_d.vertex_curves(v1).begin(), in_d.vertex_curves(v1).end());
    out_d.add_vertex_curves(new_v, other_d.vertex_curves(v2).begin(), other_d.vertex_curves(v2).end());
    return;
  }

  const Point_2& p_v = other_d.point(v);
  const X_monotone_curve_2& cv_e = in_d.edge_curves(e).front();

  Comparison_result res = m_traits->compare_y_at_x_2_object()(p_v, cv_e);

  if ((res == EQUAL) || (m_env_type == LOWER && res == SMALLER) || (m_env_type == UPPER && res == LARGER)) {
    Vertex_handle new_v = _append_vertex(out_d, p_v, e, in_d);
    out_d.add_vertex_curves(new_v, other_d.vertex_curves(v).begin(), other_d.vertex_curves(v).end());
    if (res == EQUAL) out_d.add_vertex_curves(new_v, in_d.edge_curves(e).begin(), in_d.edge_curves(e).end());
  }
}

// ---------------------------------------------------------------------------
// Compare y-coordinates of two curves at endpoints.
//
template <typename Traits, typename Diagram>
Comparison_result Envelope_divide_and_conquer_2<Traits,Diagram>::
compare_y_at_end(const X_monotone_curve_2& xcv1, const X_monotone_curve_2& xcv2, Arr_curve_end curve_end) const {
  CGAL_precondition(traits->is_in_x_range_2_object()(xcv1, xcv2));

  const auto compare_y_at_x = m_traits->compare_y_at_x_2_object();
  const auto min_vertex = m_traits->construct_min_vertex_2_object();
  const auto max_vertex = m_traits->construct_max_vertex_2_object();
  const auto param_space_in_x = m_traits->parameter_space_in_x_2_object();

  const Arr_parameter_space ps_x1 = param_space_in_x(xcv1, curve_end);
  const Arr_parameter_space ps_x2 = param_space_in_x(xcv2, curve_end);

  if (ps_x1 != ARR_INTERIOR) {
    if (ps_x2 != ARR_INTERIOR) {
      const auto cmp_y_near_boundary = m_traits->compare_y_near_boundary_2_object();
      return cmp_y_near_boundary(xcv1, xcv2, curve_end);
    }

    const auto param_space_in_y = m_traits->parameter_space_in_y_2_object();
    const Arr_parameter_space ps_y2 = param_space_in_y(xcv2, curve_end);

    if (ps_y2 == ARR_BOTTOM_BOUNDARY) return LARGER;
    else if (ps_y2 == ARR_TOP_BOUNDARY) return SMALLER;

    Comparison_result res = (curve_end == ARR_MIN_END) ?
      compare_y_at_x(min_vertex(xcv2), xcv1) :
      compare_y_at_x(max_vertex(xcv2), xcv1);

    return CGAL::opposite(res);
  }
  else if (ps_x2 != ARR_INTERIOR) {
    const auto param_space_in_y = m_traits->parameter_space_in_y_2_object();
    const Arr_parameter_space ps_y1 = param_space_in_y(xcv1, curve_end);

    if (ps_y1 == ARR_BOTTOM_BOUNDARY) return SMALLER;
    else if (ps_y1 == ARR_TOP_BOUNDARY) return LARGER;

    Comparison_result res = (curve_end == ARR_MIN_END) ?
      compare_y_at_x(min_vertex(xcv1), xcv2) :
      compare_y_at_x(max_vertex(xcv1), xcv2);
    return res;
  }

  const auto param_space_in_y = m_traits->parameter_space_in_y_2_object();
  const Arr_parameter_space ps_y1 = param_space_in_y(xcv1, curve_end);
  const Arr_parameter_space ps_y2 = param_space_in_y(xcv2, curve_end);

  if (ps_y1 != ARR_INTERIOR) {
    if (ps_y2 != ARR_INTERIOR) {
      if ((ps_y1 == ARR_BOTTOM_BOUNDARY) && (ps_y2 == ARR_TOP_BOUNDARY)) return SMALLER;
      else if ((ps_y1 == ARR_TOP_BOUNDARY) && (ps_y2 == ARR_BOTTOM_BOUNDARY)) return LARGER;

      const auto cmp_x_curve_ends = m_traits->compare_x_curve_ends_2_object();
      Comparison_result l_res = cmp_x_curve_ends(xcv1, curve_end, xcv2, curve_end);
      CGAL_assertion(l_res != EQUAL);

      if (ps_y1 == ARR_TOP_BOUNDARY) return l_res;
      else return CGAL::opposite(l_res);
    }

    const Point_2& left2 = (curve_end == ARR_MIN_END) ? min_vertex(xcv2) : max_vertex(xcv2);
    const auto cmp_x_point_curve_end = m_traits->compare_x_point_curve_end_2_object();
    Comparison_result l_res = cmp_x_point_curve_end(left2, xcv1, curve_end);

    if (l_res == LARGER) {
      Comparison_result res = compare_y_at_x(left2, xcv1);
      return CGAL::opposite(res);
    }
    else return ((ps_y1 == ARR_BOTTOM_BOUNDARY) ? SMALLER : LARGER);
  }
  else if (ps_y2 != ARR_INTERIOR) {
    const Point_2& left1 = (curve_end == ARR_MIN_END) ? min_vertex(xcv1) : max_vertex(xcv1);
    const auto cmp_x_point_curve_end = m_traits->compare_x_point_curve_end_2_object();
    Comparison_result l_res = cmp_x_point_curve_end(left1, xcv2, curve_end);
    return ((l_res == LARGER) ? compare_y_at_x(left1, xcv2) : ((ps_y2 == ARR_BOTTOM_BOUNDARY) ? LARGER : SMALLER));
  }

  const auto compare_xy = m_traits->compare_xy_2_object();
  const Point_2& left1 = (curve_end == ARR_MIN_END) ? min_vertex(xcv1) : max_vertex(xcv1);
  const Point_2& left2 = (curve_end == ARR_MIN_END) ? min_vertex(xcv2) : max_vertex(xcv2);
  Comparison_result l_res = compare_xy(left1, left2);
  return ((l_res != SMALLER) ? compare_y_at_x(left1, xcv2) : CGAL::opposite(compare_y_at_x(left2, xcv1)));
}

// ---------------------------------------------------------------------------
// Merge two non-empty intervals into the merged diagram.
//
template <typename Traits, typename Diagram>
void Envelope_divide_and_conquer_2<Traits,Diagram>::
_merge_two_intervals(Edge_const_handle e1, bool is_leftmost1, Edge_const_handle e2, bool is_leftmost2,
                     Vertex_const_handle v, bool v_exists, Comparison_result origin_of_v,
                     const Envelope_diagram_1& d1, const Envelope_diagram_1& d2, Envelope_diagram_1& out_d) {
  using Intersection_point = std::pair<Point_2, typename Traits::Multiplicity>;

  Comparison_result current_res = compare_y_at_end(d1.edge_curves(e1).front(), d2.edge_curves(e2).front(), ARR_MIN_END);
  if (m_env_type == UPPER) current_res = CGAL::opposite(current_res);

  std::optional<Point_2> p_leftmost;

  if (is_leftmost1 == true) {
    if (is_leftmost2 == false) p_leftmost = d2.point(d2.left_vertex(e2));
  }
  else {
    if (is_leftmost2 == true) p_leftmost = d1.point(d1.left_vertex(e1));
    else {
      Point_2 p1 = d1.point(d1.left_vertex(e1));
      Point_2 p2 = d2.point(d2.left_vertex(e2));
      if (m_traits->compare_xy_2_object()(p1, p2) == LARGER) p_leftmost = p1;
      else p_leftmost = p2;
    }
  }

  std::list<std::variant<Intersection_point, X_monotone_curve_2>> objects;
  const X_monotone_curve_2* intersection_curve;
  const Intersection_point* intersection_point;

  m_traits->intersect_2_object()(d1.edge_curves(e1).front(), d2.edge_curves(e2).front(), std::back_inserter(objects));

  bool equal_at_v = false;

  const Point_2* p_v = v_exists ? &( (origin_of_v == SMALLER) ? d1.point(v) : d2.point(v) ) : nullptr;

  while (! objects.empty()) {
    auto obj = std::move(objects.front());
    objects.pop_front();

    if ((intersection_point = std::get_if<Intersection_point>(&obj)) != nullptr) {
      bool is_in_x_range = true;
      if (p_leftmost && m_traits->compare_xy_2_object() (intersection_point->first, *p_leftmost) != LARGER)
        is_in_x_range = false;

      if (is_in_x_range && v_exists) {
        Comparison_result res = m_traits->compare_xy_2_object()(intersection_point->first, *p_v);

        if (res == EQUAL) equal_at_v = true;

        if (res == LARGER) break;
      }

      if (is_in_x_range) {
        CGAL_assertion(current_res != EQUAL);

        Vertex_handle new_v = (current_res == SMALLER) ?
          _append_vertex(out_d, intersection_point->first, e1, d1) :
          _append_vertex(out_d, intersection_point->first, e2, d2);

        if (equal_at_v == false) {
          out_d.add_vertex_curves(new_v, d1.edge_curves(e1).begin(), d1.edge_curves(e1).end());
          out_d.add_vertex_curves(new_v, d2.edge_curves(e2).begin(), d2.edge_curves(e2).end());
        }

        p_leftmost = intersection_point->first;
      }

      if ((current_res != EQUAL) && (intersection_point->second % 2 == 1)) {
        current_res = CGAL::opposite(current_res);
      }
      else if (((intersection_point->second == 0) || (current_res == EQUAL)) && (equal_at_v == false)) {
        current_res = m_traits->compare_y_at_x_right_2_object()(d1.edge_curves(e1).front(), d2.edge_curves(e2).front(),
                                                                intersection_point->first);
        if (m_env_type == UPPER) current_res = CGAL::opposite (current_res);
      }
    }
    else {
      intersection_curve = std::get_if<X_monotone_curve_2>(&obj);

      if (intersection_curve == nullptr) CGAL_error_msg("unrecognized intersection object.");

      auto parameter_space_x = m_traits->parameter_space_in_x_2_object();
      const bool has_left = (parameter_space_x(*intersection_curve, ARR_MIN_END) == ARR_INTERIOR);
      const bool has_right = (parameter_space_x(*intersection_curve, ARR_MAX_END) == ARR_INTERIOR);

      Point_2 pt_left, pt_right;
      if (has_left) pt_left = m_traits->construct_min_vertex_2_object()(*intersection_curve);
      if (has_right) pt_right = m_traits->construct_max_vertex_2_object()(*intersection_curve);

      bool is_in_x_range = true;
      if (p_leftmost && has_right && (m_traits->compare_xy_2_object()(pt_right, *p_leftmost) != LARGER))
        is_in_x_range = false;

      if (is_in_x_range && v_exists && has_left) {
        Comparison_result res = m_traits->compare_xy_2_object()(pt_left, *p_v);

        if (res == EQUAL) equal_at_v = true;

        if (res == LARGER) break;
      }

      if (is_in_x_range && has_left &&
          (! p_leftmost || (m_traits->compare_xy_2_object()(pt_left, *p_leftmost) == LARGER))) {
        CGAL_assertion(current_res != EQUAL);

        Vertex_handle new_v = (current_res == SMALLER) ?
          _append_vertex(out_d, pt_left, e1, d1) : _append_vertex(out_d, pt_left, e2, d2);

        if (equal_at_v == false) {
          out_d.add_vertex_curves(new_v, d1.edge_curves(e1).begin(), d1.edge_curves(e1).end());
          out_d.add_vertex_curves(new_v, d2.edge_curves(e2).begin(), d2.edge_curves(e2).end());
        }

        p_leftmost = pt_left;
      }

      if (is_in_x_range && has_right && (! v_exists || (m_traits->compare_xy_2_object()(pt_right, *p_v) == SMALLER))) {
        Vertex_handle new_v = _append_vertex(out_d, pt_right, e1, d1);
        Edge_handle e_left = out_d.left_edge(new_v);
        out_d.add_edge_curves(e_left, d2.edge_curves(e2).begin(), d2.edge_curves(e2).end());

        CGAL_assertion(equal_at_v == false);
        out_d.add_vertex_curves(new_v, d1.edge_curves(e1).begin(), d1.edge_curves(e1).end());
        out_d.add_vertex_curves(new_v, d2.edge_curves(e2).begin(), d2.edge_curves(e2).end());

        p_leftmost = pt_right;
      }

      if (has_right == false || (v_exists && m_traits->compare_xy_2_object()(pt_right, *p_v) != SMALLER)) {
        if (v_exists) {
          Vertex_handle new_v = _append_vertex(out_d, *p_v, e1, d1);
          Edge_handle e_left = out_d.left_edge(new_v);
          out_d.add_edge_curves(e_left, d2.edge_curves(e2).begin(), d2.edge_curves(e2).end());
        }

        equal_at_v = true;
        current_res = EQUAL;
        break;
      }

      current_res =
        m_traits->compare_y_at_x_right_2_object()(d1.edge_curves(e1).front(), d2.edge_curves(e2).front(), pt_right);
      if (m_env_type == UPPER) current_res = CGAL::opposite (current_res);
    }
  }

  if (equal_at_v) {
    CGAL_assertion(v_exists);

    Edge_handle e_last = out_d.rightmost();
    if (out_d.has_left_vertex(e_last)) {
      Vertex_handle v_to_be_updated = out_d.left_vertex(e_last);
      if (origin_of_v == EQUAL) {
        Vertex_const_handle v1 = d1.right_vertex(e1);
        Vertex_const_handle v2 = d2.right_vertex(e2);
        out_d.add_vertex_curves(v_to_be_updated, d1.vertex_curves(v1).begin(), d1.vertex_curves(v1).end());
        out_d.add_vertex_curves(v_to_be_updated, d2.vertex_curves(v2).begin(), d2.vertex_curves(v2).end());
      }
      else {
        const Envelope_diagram_1& src_d = (origin_of_v == SMALLER) ? d2 : d1;
        Edge_const_handle e = (origin_of_v == SMALLER) ? e2 : e1;
        out_d.add_vertex_curves(v_to_be_updated, src_d.vertex_curves(v).begin(), src_d.vertex_curves(v).end());
        out_d.add_vertex_curves(v_to_be_updated, src_d.edge_curves(e).begin(), src_d.edge_curves(e).end());
      }
    }
    return;
  }

  if (! v_exists) {
    switch (current_res) {
     case SMALLER:
      out_d.add_edge_curves(out_d.rightmost(), d1.edge_curves(e1).begin(), d1.edge_curves(e1).end());
      return;

     case LARGER:
      out_d.add_edge_curves(out_d.rightmost(), d2.edge_curves(e2).begin(), d2.edge_curves(e2).end());
      return;

     case EQUAL:
      out_d.add_edge_curves(out_d.rightmost(), d1.edge_curves(e1).begin(), d1.edge_curves(e1).end());
      out_d.add_edge_curves(out_d.rightmost(), d2.edge_curves(e2).begin(), d2.edge_curves(e2).end());
      return;

     default:
      CGAL_error_msg("should not reach here.");
      return;
    }
  }

  if (current_res == SMALLER) {
    if (origin_of_v == SMALLER) {
      Vertex_handle new_v = _append_vertex(out_d, d1.point(v), e1, d1);
      out_d.add_vertex_curves(new_v, d1.vertex_curves(v).begin(), d1.vertex_curves(v).end());
      return;
    }

    if (origin_of_v == EQUAL) {
      Vertex_handle new_v = _append_vertex(out_d, d2.point(v), e1, d1);
      out_d.add_vertex_curves(new_v, d2.vertex_curves(v).begin(), d2.vertex_curves(v).end());
      Vertex_const_handle v1 = d1.right_vertex(e1);
      out_d.add_vertex_curves(new_v, d1.vertex_curves(v1).begin(), d1.vertex_curves(v1).end());
      return;
    }

    const Point_2& p2 = d2.point(v);
    const Comparison_result res = m_traits->compare_y_at_x_2_object()(p2, d1.edge_curves(e1).front());
    if (res == EQUAL || ((m_env_type == LOWER) && (res == SMALLER)) || ((m_env_type == UPPER) && (res == LARGER))) {
      Vertex_handle new_v = _append_vertex(out_d, p2, e1, d1);
      out_d.add_vertex_curves(new_v, d2.vertex_curves(v).begin(), d2.vertex_curves(v).end());
      if (res == EQUAL) out_d.add_vertex_curves(new_v, d1.edge_curves(e1).begin(), d1.edge_curves(e1).end());
    }
    return;
  }

  if (origin_of_v != SMALLER) {
    Vertex_handle new_v = _append_vertex(out_d, d2.point(v), e2, d2);
    out_d.add_vertex_curves(new_v, d2.vertex_curves(v).begin(), d2.vertex_curves(v).end());
    if (origin_of_v == EQUAL) {
      Vertex_const_handle v1 = d1.right_vertex(e1);
      out_d.add_vertex_curves(new_v, d1.vertex_curves(v1).begin(), d1.vertex_curves(v1).end());
    }
    return;
  }

  const Point_2& p1 = d1.point(v);
  const Comparison_result res = m_traits->compare_y_at_x_2_object()(p1, d2.edge_curves(e2).front());
  if (res == EQUAL || ((m_env_type == LOWER) && (res == SMALLER)) || ((m_env_type == UPPER) && (res == LARGER))) {
    Vertex_handle new_v = _append_vertex(out_d, p1, e2, d2);
    out_d.add_vertex_curves(new_v, d1.vertex_curves(v).begin(), d1.vertex_curves(v).end());
    if (res == EQUAL) out_d.add_vertex_curves(new_v, d2.edge_curves(e2).begin(), d2.edge_curves(e2).end());
  }
}

// ---------------------------------------------------------------------------
// Append a vertex to the given diagram.
//
template <typename Traits, typename Diagram>
typename Envelope_divide_and_conquer_2<Traits,Diagram>::Vertex_handle
Envelope_divide_and_conquer_2<Traits,Diagram>::
_append_vertex(Envelope_diagram_1& diag, const Point_2& p, Edge_const_handle e, const Envelope_diagram_1& in_d) {
  auto& topo = diag.topology_traits();

  // Create the new vertex and the new edge.
  Vertex_handle new_v = diag.new_vertex(p);
  Edge_handle new_e = diag.new_edge();

  if (! in_d.empty_edge_curves(e)) diag.add_edge_curves(new_e, in_d.edge_curves(e).begin(), in_d.edge_curves(e).end());

  // Connect the new vertex.
  topo.set_left_edge(new_v, new_e);
  topo.set_right_edge(new_v, diag.rightmost());

  if (diag.leftmost() != diag.rightmost()) {
    // The diagram is not empty. Connect the new edge to the left of the
    // rightmost edge of the diagram.
    topo.set_right_vertex(new_e, new_v);

    Vertex_handle last_v = diag.left_vertex(diag.rightmost());
    topo.set_left_vertex(new_e, last_v);
    topo.set_right_edge(last_v, new_e);
    topo.set_left_vertex(diag.rightmost(), new_v);
  }
  else {
    // The diagram is empty: Make the new edge the leftmost.
    topo.set_right_vertex(new_e, new_v);
    diag.set_leftmost(new_e);
    topo.set_left_vertex(diag.rightmost(), new_v);
  }

  return new_v;
}

// ---------------------------------------------------------------------------
// Merge the vertical segments into the envelope given as a minimization
// (or maximization) diagram.
//
template <typename Traits, typename Diagram>
void Envelope_divide_and_conquer_2<Traits, Diagram>::
_merge_vertical_segments(Curve_pointer_vector& vert_vec, Envelope_diagram_1& out_d) {
  Less_vertical_segment les_vert(m_traits);

  std::sort(vert_vec.begin(), vert_vec.end(), les_vert);

  typename Traits_adaptor_2::Compare_x_2 comp_x = m_traits->compare_x_2_object();
  typename Traits_adaptor_2::Compare_xy_2 comp_xy = m_traits->compare_xy_2_object();
  typename Traits_adaptor_2::Compare_y_at_x_2 comp_y_at_x = m_traits->compare_y_at_x_2_object();
  typename Traits_adaptor_2::Construct_min_vertex_2 min_vertex = m_traits->construct_min_vertex_2_object();
  typename Traits_adaptor_2::Construct_max_vertex_2 max_vertex = m_traits->construct_max_vertex_2_object();

  Edge_handle e = out_d.leftmost();
  Vertex_handle v{};
  Curve_pointer_iterator iter = vert_vec.begin();
  Curve_pointer_iterator next;
  Comparison_result res;
  bool in_e_range;
  bool on_v;
  Point_2 p;

  while (iter != vert_vec.end()) {
    if (e != out_d.rightmost()) {
      v = out_d.right_vertex(e);

      res = comp_x(min_vertex(**iter), out_d.point(v));
      in_e_range = (res != LARGER);
      on_v = (res == EQUAL);
    }
    else {
      in_e_range = true;
      on_v = false;
    }

    if (! in_e_range) {
      e = out_d.right_edge(v);
      continue;
    }

    std::list<X_monotone_curve_2> env_cvs;

    env_cvs.push_back(**iter);
    next = iter;
    ++next;
    while ((next != vert_vec.end()) && (comp_x(min_vertex(**iter), min_vertex(**next)) == EQUAL)) {
      if (m_env_type == LOWER) {
        res = comp_xy(min_vertex(env_cvs.front()), min_vertex(**next));

        if (res == EQUAL) {
          env_cvs.push_back(**next);
        }
        if (res == LARGER) {
          env_cvs.clear();
          env_cvs.push_back(**next);
        }
      }
      else {
        res = comp_xy(max_vertex(env_cvs.front()), max_vertex(**next));

        if (res == EQUAL) {
          env_cvs.push_back(**next);
        }
        if (res == SMALLER) {
          env_cvs.clear();
          env_cvs.push_back(**next);
        }
      }

      ++next;
    }

    p = (m_env_type == LOWER) ?
      min_vertex(env_cvs.front()) : max_vertex(env_cvs.front());

    if (on_v) {
      res = comp_xy(p, out_d.point(v));

      if (res == EQUAL) {
        out_d.add_vertex_curves(v, env_cvs.begin(), env_cvs.end());
      }
      else if ((m_env_type == LOWER && res == SMALLER) || (m_env_type == UPPER && res == LARGER)) {
        out_d.clear_vertex_curves(v);
        out_d.add_vertex_curves(v, env_cvs.begin(), env_cvs.end());
      }
    }
    else {
      if (out_d.empty_edge_curves(e)) {
        Vertex_handle new_v = _split_edge(out_d, p, e);
        out_d.add_vertex_curves(new_v, env_cvs.begin(), env_cvs.end());
      }
      else {
        res = comp_y_at_x(p, out_d.edge_curves(e).front());
        if (((m_env_type == LOWER) && (res != LARGER)) || ((m_env_type == UPPER) && (res != SMALLER))) {
          Vertex_handle new_v = _split_edge(out_d, p, e);
          out_d.add_vertex_curves(new_v, env_cvs.begin(), env_cvs.end());
          if (res == EQUAL) out_d.add_vertex_curve(new_v, out_d.edge_curves(e).front());
        }
      }
    }

    iter = next;
  }
}

// ---------------------------------------------------------------------------
// Split a given diagram edge by inserting a vertex in its interior.
//
template <typename Traits, typename Diagram>
typename Envelope_divide_and_conquer_2<Traits,Diagram>::Vertex_handle
Envelope_divide_and_conquer_2<Traits,Diagram>::
_split_edge(Envelope_diagram_1& diag, const Point_2& p, Edge_handle e) {
  auto& topo = diag.topology_traits();

  // Create the new vertex and the new edge.
  Vertex_handle new_v = diag.new_vertex(p);
  Edge_handle new_e = diag.new_edge();

  // Duplicate the curves container associated with e.
  if (! diag.empty_edge_curves(e)) diag.add_edge_curves(new_e, diag.edge_curves(e).begin(), diag.edge_curves(e).end());

  // Connect the new vertex between e and new_e.
  topo.set_left_edge(new_v, e);
  topo.set_right_edge(new_v, new_e);

  topo.set_left_vertex(new_e, new_v);
  if (e != diag.rightmost()) {
    Vertex_handle v_right = diag.right_vertex(e);
    topo.set_right_vertex(new_e, v_right);
    topo.set_left_edge(v_right, new_e);
  }
  else {
    diag.set_rightmost(new_e);
  }

  topo.set_right_vertex(e, new_v);

  return new_v;
}

} // namespace CGAL

#endif
