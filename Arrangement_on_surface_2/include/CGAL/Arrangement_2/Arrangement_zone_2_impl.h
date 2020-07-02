// Copyright (c) 2005,2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s): Ron Wein          <wein@post.tau.ac.il>
//            Efi Fogel         <efif@post.tau.ac.il>
//            (based on old version by Eyal Flato)

#ifndef CGAL_ARRANGEMENT_ZONE_2_IMPL_H
#define CGAL_ARRANGEMENT_ZONE_2_IMPL_H

#include <CGAL/license/Arrangement_on_surface_2.h>

/*! \file
 * Member-function definitions for the Arrangement_zone_2 class.
 */

namespace CGAL {

//-----------------------------------------------------------------------------
// Initialize the zone-computation process with a given curve and an object
// that wraps the location of the curve's left end.
//
template <typename Arrangement, typename ZoneVisitor>
void Arrangement_zone_2<Arrangement, ZoneVisitor>::
init_with_hint(const X_monotone_curve_2& cv, Pl_result_type obj)
{
  // Set the curve and check whether its ends are bounded, therefore
  // associated with valid endpoints.
  m_cv = cv;

  if (m_geom_traits->is_closed_2_object()(m_cv, ARR_MIN_END)) {
    // The left endpoint is valid.
    const Arr_parameter_space ps_x1 =
      m_geom_traits->parameter_space_in_x_2_object()(m_cv, ARR_MIN_END);
    const Arr_parameter_space ps_y1 =
      m_geom_traits->parameter_space_in_y_2_object()(m_cv, ARR_MIN_END);
    m_has_left_pt = true;
    m_left_on_boundary = (ps_x1 != ARR_INTERIOR || ps_y1 != ARR_INTERIOR);
    m_left_pt = m_geom_traits->construct_min_vertex_2_object()(m_cv);
  }
  else {
    // The left end of the curve lies on open boundary.
    m_has_left_pt = false;
    m_left_on_boundary = true;
  }

  if (m_geom_traits->is_closed_2_object()(m_cv, ARR_MAX_END)) {
    // The right endpoint is valid.
    const Arr_parameter_space ps_x2 =
      m_geom_traits->parameter_space_in_x_2_object()(m_cv, ARR_MAX_END);
    const Arr_parameter_space ps_y2 =
      m_geom_traits->parameter_space_in_y_2_object()(m_cv, ARR_MAX_END);
    m_has_right_pt = true;
    m_right_on_boundary = (ps_x2 != ARR_INTERIOR || ps_y2 != ARR_INTERIOR);
    m_right_pt = m_geom_traits->construct_max_vertex_2_object()(m_cv);
  }
  else {
    // The right end of the curve lies on open boundary.
    m_has_right_pt = false;
    m_right_on_boundary = true;
  }

  // Set the object that represents the location of the left end of the curve
  // in the arrangement.
  m_obj = obj;
}

//-----------------------------------------------------------------------------
// Compute the zone of the given curve and issue the apporpriate
// notifications for the visitor.
//
template <typename Arrangement, typename ZoneVisitor>
void Arrangement_zone_2<Arrangement, ZoneVisitor>::compute_zone()
{
  // Initialize flags and set all handles to be invalid.
  bool done = false;

  m_found_intersect = false;
  m_found_overlap = false;
  m_found_iso_vert = false;

  m_left_v = m_invalid_v;
  m_left_he = m_invalid_he;
  m_right_v = m_invalid_v;
  m_right_he = m_invalid_he;

  // Locate the arrangement feature containing the left endpoint of the
  // curve (currently m_obj stores the object containing it).
  const Vertex_const_handle* vh;
  const Halfedge_const_handle* hh;

  if ((vh = boost::get<Vertex_const_handle>(&m_obj)) != nullptr) {
    CGAL_assertion(m_has_left_pt);

    // The left endpoint coincides with an existing vertex:
    m_left_v = m_arr.non_const_handle(*vh);

#if 0
    if (m_left_on_boundary) {
      // Use the accessor to locate the predecessor edge, in case the left
      // endpoint has boundary conditions.
      const Arr_parameter_space ps_x =
        m_geom_traits->parameter_space_in_x_2_object()(m_cv, ARR_MIN_END);
      const Arr_parameter_space ps_y =
        m_geom_traits->parameter_space_in_y_2_object()(m_cv, ARR_MIN_END);

      m_left_he = m_arr_access.locate_around_boundary_vertex(m_left_v,
                                                             m_cv, ARR_MIN_END,
                                                             ps_x, ps_y);
    }
#endif

  }
  else if ((hh = boost::get<Halfedge_const_handle>(&m_obj)) != nullptr) {
    if (m_has_left_pt) {
      // Obtain the right halfedge from the halfedge-pair containing m_left_pt
      // in their interior.
      m_left_he =
        _direct_intersecting_edge_to_right(m_cv, m_left_pt,
                                           m_arr.non_const_handle(*hh));

      // Handle overlaps.
      if (m_found_overlap) {
        // In this case m_cv overlaps the curve associated with m_intersect_he.
        // Compute the overlapping subcurve.
        bool dummy;
        auto obj = _compute_next_intersection(m_intersect_he, false, dummy);
        m_overlap_cv = boost::get<X_monotone_curve_2>(*obj);

        // Remove the overlap from the map.
        _remove_next_intersection(m_intersect_he);

        // Compute the overlap zone.
        done = _zone_in_overlap();
      }
    }
    else {
      // In case the unbounded left end conicides with an edge, then our curve
      // overlaps the curve associated with this edge.
      m_intersect_he = m_arr.non_const_handle(*hh);

      bool dummy;
      auto obj = _compute_next_intersection(m_intersect_he, false, dummy);
      m_overlap_cv = boost::get<X_monotone_curve_2>(*obj);

      // Remove the overlap from the map.
      _remove_next_intersection(m_intersect_he);

      // Compute the overlap zone.
      done = _zone_in_overlap();
    }
  }
  else {
    // The left endpoint lies inside a face.
    const Face_const_handle* fh = boost::get<Face_const_handle>(&m_obj);

    CGAL_assertion_msg(fh != nullptr,
                       "Invalid object returned by the point-location query.");

    // Compute the zone of the curve at the interior of the face.
    // m_left_pt is not on the face boundary.
    done = _zone_in_face(m_arr.non_const_handle(*fh), false);

    // In case we have just discovered an overlap, compute the overlapping
    // zone as well.
    if (! done && m_found_overlap) done = _zone_in_overlap();
  }

  // Compute the zone of the curve (or what is remaining of it) in the
  // arrangement, starting from the current position we have computed.
  while (! done) {
    // Check if we know the face the curve is going to penetrate now.
    if (m_left_he == m_invalid_he) {
      if (m_left_v != m_invalid_v) {
        // We know the vertex that coincides with the left endpoint of m_cv.
        if (! m_left_v->is_isolated()) {
          // Locate the curve around the m_left_v vertex - that is, find a
          // halfedge m_left_he such that m_cv should be placed between
          // m_left_he and its current successor around the vertex, going in a
          // clockwise order.
          m_found_overlap = _find_prev_around_vertex(m_left_v, m_left_he);
        }
        else {
          // m_left_v is an isolated vertex.
          m_found_iso_vert = true;
        }
      }
      else {
        CGAL_assertion(m_right_he != m_invalid_he);

        // In this case m_right_he is the halfedge that the left portion of m_cv
        // intersected, and we obtain m_left_he by comparing the remaining
        // portion of m_cv with the curve associated with this edge.
        m_left_he =
          _direct_intersecting_edge_to_right(m_cv, m_left_pt, m_right_he);
      }

      if (m_found_overlap) {
        // In this case m_cv overlaps the curve associated with m_intersect_he.
        // Compute the overlapping subcurve to the right of curr_v.
        bool dummy;
        auto obj = _compute_next_intersection(m_intersect_he, false, dummy);
        m_overlap_cv = boost::get<X_monotone_curve_2>(*obj);

        // Remove the overlap from the map.
        _remove_next_intersection(m_intersect_he);

        // Compute the overlap zone and continue to the end of the loop.
        done = _zone_in_overlap();
        continue;
      }
    }

    if ((m_left_v == m_invalid_v) || ! m_left_v->is_isolated()) {
      // At this point we can compute the zone of m_cv starting from the
      // m_left_he inside its incident face.
      done = _zone_in_face(m_left_he->face(), true);
      // m_left_pt is on the face boundary.
    }
    else {
      // Compute the zone of m_cv starting from the face that contains the
      // isolated vertex m_left_v.
      done = _zone_in_face(m_left_v->face(), false);
      // m_left_pt is not on the face boundary.
    }

    // In case we have just discovered an overlap, compute the overlapping
    // zone as well.
    if (! done && m_found_overlap) done = _zone_in_overlap();
  }

  // Clear the intersections map.
  m_inter_map.clear();
}

//-----------------------------------------------------------------------------
// Check whether two curves with a common min endpoint overlap.
//
template <typename Arrangement, typename ZoneVisitor>
bool  Arrangement_zone_2<Arrangement, ZoneVisitor>::
do_overlap_impl(const X_monotone_curve_2& cv1,
                const X_monotone_curve_2& cv2,
                const Point_2& p, Arr_not_all_sides_oblivious_tag) const
{
  typename Traits_adaptor_2::Compare_y_at_x_right_2 cmp_right =
    m_geom_traits->compare_y_at_x_right_2_object();

  Arr_parameter_space psx =
    m_geom_traits->parameter_space_in_x_2_object()(cv1, ARR_MIN_END);
  Arr_parameter_space psy =
    m_geom_traits->parameter_space_in_y_2_object()(cv1, ARR_MIN_END);

  if ((psx == ARR_INTERIOR) && (psy == ARR_INTERIOR))
    return (cmp_right(cv1, cv2, p) == EQUAL);

  bool vertical1 = m_geom_traits->is_vertical_2_object()(cv1);
  bool vertical2 = m_geom_traits->is_vertical_2_object()(cv2);
  if (vertical1 != vertical2) return false;

  if (psx == ARR_LEFT_BOUNDARY) {
    // If both curves are vertical and their common min endpoint lies on the
    // left boundary, they completely lie on the left boundary and they overlap.
    if (vertical1) return true;

    typename Traits_adaptor_2::Compare_y_near_boundary_2 cmp_near =
      m_geom_traits->compare_y_near_boundary_2_object();
    return (cmp_near(cv1, cv2, ARR_MIN_END) == EQUAL);
  }

  // If the common min endpoint lies on the right boundary, both curves are
  // vertical, they completely lie on the right boundary, and they overlap.
  if (psx == ARR_RIGHT_BOUNDARY) return true;

  // If the curves are not vertical, we can safly call the standard function.
  // Observe that this case covers the case where (psy == ARR_TOP_BOUNDARY).
  if (! vertical1) return (cmp_right(cv1, cv2, p) == EQUAL);

  CGAL_assertion(psy == ARR_BOTTOM_BOUNDARY);
  // Observe, that if both curves are vertical and have a common min endpoint
  // that lies on the bottom boundary, and the bottom boundary is contracted
  // the curves do not necessarily overlap.
  // In this case it is sufficient to test whether the x-coordinates on the
  // boundary are equal.
  typename Traits_adaptor_2::Compare_x_on_boundary_2 cmp_x_on_bd =
    m_geom_traits->compare_x_on_boundary_2_object();
  typename Traits_adaptor_2::Compare_x_near_boundary_2 cmp_x_near_bd =
    m_geom_traits->compare_x_near_boundary_2_object();
  Comparison_result res = cmp_x_on_bd(cv1, ARR_MIN_END, cv2, ARR_MIN_END);
  if (res == EQUAL) res = cmp_x_near_bd(cv1, cv2, ARR_MIN_END);
  return (res == EQUAL);
}

//-----------------------------------------------------------------------------
// Check whether the given query curve is encountered when rotating the
// first curve in a clockwise direction around a given point until reaching
// the second curve.
template <typename Arrangement, typename ZoneVisitor>
bool Arrangement_zone_2<Arrangement, ZoneVisitor>::
is_between_cw_impl(const X_monotone_curve_2& cv, bool cv_to_right,
                   const X_monotone_curve_2& cv1, bool cv1_to_right,
                   const X_monotone_curve_2& cv2, bool cv2_to_right,
                   const Point_2& p,
                   bool& cv_equal_cv1, bool& cv_equal_cv2,
                   Arr_not_all_sides_oblivious_tag) const
{
  cv_equal_cv1 = false;
  cv_equal_cv2 = false;

  typename Traits_adaptor_2::Parameter_space_in_x_2 ps_in_x =
    m_geom_traits->parameter_space_in_x_2_object();
  typename Traits_adaptor_2::Parameter_space_in_y_2 ps_in_y =
    m_geom_traits->parameter_space_in_y_2_object();
  typename Traits_adaptor_2::Is_between_cw_2 is_bcw =
    m_geom_traits->is_between_cw_2_object();

  bool vertical = m_geom_traits->is_vertical_2_object()(cv);
  bool vertical1 = m_geom_traits->is_vertical_2_object()(cv1);
  bool vertical2 = m_geom_traits->is_vertical_2_object()(cv2);

  typename Traits_adaptor_2::Compare_x_on_boundary_2 cmp_x_on_bd =
    m_geom_traits->compare_x_on_boundary_2_object();
  typename Traits_adaptor_2::Compare_x_near_boundary_2 cmp_x_near_bd =
    m_geom_traits->compare_x_near_boundary_2_object();
  typename Traits_adaptor_2::Compare_y_near_boundary_2 cmp_y_near_bd =
    m_geom_traits->compare_y_near_boundary_2_object();

  ////
  if (!cv1_to_right) {
    Arr_parameter_space psx = ps_in_x(cv1, ARR_MAX_END);
    Arr_parameter_space psy = ps_in_y(cv1, ARR_MAX_END);
    if ((psx == ARR_INTERIOR) && (psy == ARR_INTERIOR))
      return is_bcw(cv, cv_to_right, cv1, cv1_to_right, cv2, cv2_to_right,
                    p, cv_equal_cv1, cv_equal_cv2);

    if (!cv2_to_right) {
      // Case 1: Both cv1 and cv2 are defined to the left of p.

      CGAL_assertion(psx != ARR_RIGHT_BOUNDARY);

      if (psx == ARR_LEFT_BOUNDARY) {
        // Case 1.1 The point resides on the right boundary.
        CGAL_assertion(vertical);

        /* Case 1.1.1 cv1 is vertical
         *       o
         *       |cv
         *   cv2 |
         * o-----o v
         *       |
         *       |cv1
         *       o
         */
        if (vertical1) return false;

        /* Case 1.1.2 cv2 is vertical,
         *       o
         *       |cv
         *   cv1 |
         * o-----o v
         *       |
         *       |cv2
         *       o
         */
        if (vertical2) return true;

        /* Case 1.1.3
         *    o  o
         *     \ |cv
         *      \|
         *       o
         *      /
         *     /
         *    o
         */
        return (cmp_y_near_bd(cv1, cv2, ARR_MAX_END) == LARGER);
      }

      CGAL_assertion((psy == ARR_TOP_BOUNDARY) || (psy == ARR_BOTTOM_BOUNDARY));
      // Case 1.2 & 1.3 The point resides on the top or bottom boundaries
      Comparison_result res = cmp_x_on_bd(cv1, ARR_MAX_END, cv2, ARR_MAX_END);
      if (res == SMALLER) return (psy == ARR_TOP_BOUNDARY);
      if (res == LARGER) return (psy != ARR_TOP_BOUNDARY);

      res = cmp_x_on_bd(cv1, ARR_MAX_END, cv2, ARR_MAX_END);
      if (res == SMALLER) return (psy == ARR_TOP_BOUNDARY);
      CGAL_assertion(res == LARGER);
      return (psy != ARR_TOP_BOUNDARY);
    }

    CGAL_assertion(!cv1_to_right && cv2_to_right);
    /* Case 3: cv1 is defined to the left of p, and cv2 to its right.
     *   cv1   cv2
     * o-----o-----o
     */
    if (psx == ARR_LEFT_BOUNDARY) {
      /* Case 3.1
       *   cv2
       * o-----o
       * |
       * |cv1
       * o
       */
      CGAL_assertion(vertical1);
      if (vertical) return true;

      Comparison_result res = cmp_y_near_bd(cv, cv2, ARR_MIN_END);
      if (res == EQUAL) {
        cv_equal_cv2 = true;
        return false;
      }
      return (res == LARGER);
    }

    if (psx == ARR_RIGHT_BOUNDARY) {
      /* Case 3.2
       *       o
       *       |cv2
       *   cv1 |
       * o-----o
       */
      CGAL_assertion(vertical2);
      CGAL_assertion(vertical);
      cv_equal_cv2 = true;
      return false;
    }

    /* Case 3.3 & Case 4.4
     * o       o       o
     *  \     /    cv1/ \cv2
     *   \   /       /   \
     * cv1\ /cv2    /     \
     *     o       o       o
     */
    Comparison_result res2 = cmp_x_on_bd(cv, ARR_MIN_END, cv2, ARR_MIN_END);
    if (res2 == EQUAL) res2 = cmp_x_near_bd(cv, cv2, ARR_MIN_END);
    if (res2 == SMALLER) return (psy == ARR_BOTTOM_BOUNDARY);
    if (res2 == EQUAL) {
      cv_equal_cv2 = true;
      return false;
    }
    CGAL_assertion(res2 == LARGER);
    return (psy == ARR_TOP_BOUNDARY);
  }

  ////////
  Arr_parameter_space psx = ps_in_x(cv1, ARR_MIN_END);
  Arr_parameter_space psy = ps_in_y(cv1, ARR_MIN_END);
  if ((psx == ARR_INTERIOR) && (psy == ARR_INTERIOR))
    return is_bcw(cv, cv_to_right, cv1, cv1_to_right, cv2, cv2_to_right,
                  p, cv_equal_cv1, cv_equal_cv2);
  if (cv2_to_right) {
    // Case 2: Both cv1 and cv2 are defined to the right of p.
    CGAL_assertion(psx != ARR_RIGHT_BOUNDARY);
    if (psx == ARR_LEFT_BOUNDARY) {
      /* case 2.1.1 cv1 is vertical
       * o
       * |cv1
       * |
       * o------o
       *    cv2
       */
      if (vertical1) {
        if (vertical) {
          cv_equal_cv1 = true;
          return false;
        }

        Comparison_result res = cmp_y_near_bd(cv, cv2, ARR_MIN_END);
        if (res == EQUAL) {
          cv_equal_cv2 = true;
          return false;
        }
        return (res == LARGER);
      }

      /* case 2.1.2 cv2 is vertical
       * o
       * |cv2
       * |
       * o------o
       *    cv1
       */
      if (vertical2) {
        if (vertical) {
          cv_equal_cv2 = true;
          return false;
        }

        Comparison_result res = cmp_y_near_bd(cv, cv1, ARR_MIN_END);
        if (res == EQUAL) {
          cv_equal_cv1 = true;
          return false;
        }
        return (res == SMALLER);
      }

      // case 2.1.3 cv is vertical
      if (vertical) {
        Comparison_result res = cmp_y_near_bd(cv1, cv2, ARR_MIN_END);
        return (res == LARGER);
      }

      // case 2.1.4 none are vertical
      Comparison_result res1 = cmp_y_near_bd(cv, cv1, ARR_MIN_END);
      if (res1 == EQUAL) {
        cv_equal_cv1 = true;
        return false;
      }
      Comparison_result res2 = cmp_y_near_bd(cv, cv2, ARR_MIN_END);
      if (res2 == EQUAL) {
        cv_equal_cv2 = true;
        return false;
      }
      Comparison_result res = cmp_y_near_bd(cv1, cv2, ARR_MIN_END);

      if ((res1 == LARGER) && (res2 == SMALLER)) return false;
      if ((res1 == SMALLER) && (res2 == LARGER)) return true;
      if ((res1 == LARGER) && (res2 == LARGER)) return (res == SMALLER);
      CGAL_assertion((res1 == SMALLER) && (res2 == SMALLER));
      return (res == LARGER);
    }

    CGAL_assertion((psy == ARR_TOP_BOUNDARY) || (psy == ARR_BOTTOM_BOUNDARY));
    // Case 2.2 & 2.3 The point resides on the top or bottom boundaries
    Comparison_result res1 = cmp_x_on_bd(cv, ARR_MIN_END, cv1, ARR_MIN_END);
    if (res1 == EQUAL) res1 = cmp_x_near_bd(cv, cv1, ARR_MIN_END);
    Comparison_result res2 = cmp_x_on_bd(cv, ARR_MIN_END, cv2, ARR_MIN_END);
    if (res2 == EQUAL) res2 = cmp_x_near_bd(cv, cv2, ARR_MIN_END);
    Comparison_result res = cmp_x_on_bd(cv1, ARR_MIN_END, cv2, ARR_MIN_END);
    if (res == EQUAL) res = cmp_x_near_bd(cv1, cv2, ARR_MIN_END);

    if (res1 == EQUAL) {
      cv_equal_cv1 = true;
      return false;
    }
    if (res2 == EQUAL) {
      cv_equal_cv2 = true;
      return false;
    }
    if ((res1 == LARGER) && (res2 == SMALLER)) return (psy != ARR_TOP_BOUNDARY);
    if ((res1 == SMALLER) && (res2 == LARGER)) return (psy == ARR_TOP_BOUNDARY);
    if ((res1 == LARGER) && (res2 == LARGER))
      return (res == SMALLER) ?
        (psy == ARR_TOP_BOUNDARY) : (psy != ARR_TOP_BOUNDARY);
    CGAL_assertion((res1 == SMALLER) && (res2 == SMALLER));
    return (res == LARGER) ?
      (psy != ARR_TOP_BOUNDARY) : (psy == ARR_TOP_BOUNDARY);
  }

  CGAL_assertion(cv1_to_right && ! cv2_to_right);
  /* Case 4: cv1 is defined to the right of p, and cv2 to its left.
   *   cv2   cv1
   * o-----o-----o
   */
  if (psx == ARR_LEFT_BOUNDARY) {
    /* Case 4.1
     *   cv1
     * o-----o
     * |
     * |cv2
     * o
     */
    CGAL_assertion(vertical2);
    if (vertical) return true;

    Comparison_result res = cmp_y_near_bd(cv, cv1, ARR_MIN_END);
    if (res == EQUAL) {
      cv_equal_cv1 = true;
      return false;
    }
    return (res == SMALLER);
  }

  if (psx == ARR_RIGHT_BOUNDARY) {
    /* Case 4.2
     *       o
     *       |cv1
     *   cv2 |
     * o-----o
     */
    CGAL_assertion(vertical1);
    CGAL_assertion(vertical);
    cv_equal_cv1 = true;
    return false;
  }

  /* Case 4.3 & Case 4.4
   * o       o       o
   *  \     /    cv2/ \cv1
   *   \   /       /   \
   * cv2\ /cv1    /     \
   *     o       o       o
   */
  Comparison_result res1 = cmp_x_on_bd(cv, ARR_MIN_END, cv1, ARR_MIN_END);
  if (res1 == EQUAL) res1 = cmp_x_near_bd(cv, cv1, ARR_MIN_END);
  if (res1 == LARGER) return (psy == ARR_TOP_BOUNDARY);
  if (res1 == EQUAL) {
    cv_equal_cv1 = true;
    return false;
  }
  CGAL_assertion(res1 == SMALLER);
  return (psy != ARR_TOP_BOUNDARY);
}

//-----------------------------------------------------------------------------
// Find a face containing the query curve m_cv around the given vertex.
// In case an overlap occurs, sets m_intersect_he to be the overlapping edge.
//
template <typename Arrangement, typename ZoneVisitor>
bool Arrangement_zone_2<Arrangement, ZoneVisitor>::
_find_prev_around_vertex(Vertex_handle v, Halfedge_handle& he)
{

  // Go over the incident halfedges of v, going in a clockwise order.
  typename Arrangement_2::Halfedge_around_vertex_circulator he_first =
    v->incident_halfedges();
  typename Arrangement_2::Halfedge_around_vertex_circulator he_curr = he_first;
  typename Arrangement_2::Halfedge_around_vertex_circulator he_next = he_curr;
  ++he_next;

  if (he_curr == he_next) {
    // In case there is just a single incident halfedge around v,
    // we should insert m_cv right after this halfedge.
    he = he_curr;

    // Note that m_cv extends to the right of v.
    // If the single halfedge he extends to the left, there is no overlap.
    if (he->direction() != ARR_RIGHT_TO_LEFT) return false;

    // Both extends to the right; check for an overlap.
    if (! do_overlap(he->curve(), m_cv, v->point())) return false;

    // Mark that an overlap occurs:
    m_intersect_he = he_curr;
    return true;
  }

  // Find the face containing m_cv around the vertex.
  do {
    // Check if it is possible to insert cv in between the current curve
    // and the next curve, going in a clockwise direction around v.
    bool cv_equals_curr;
    bool cv_equals_next;
    bool is_between = is_between_cw(m_cv, true,
                                    he_curr->curve(),
                                    (he_curr->direction() == ARR_RIGHT_TO_LEFT),
                                    he_next->curve(),
                                    (he_next->direction() == ARR_RIGHT_TO_LEFT),
                                    v->point(),
                                    cv_equals_curr, cv_equals_next);

    // Check the case of overlaps:
    if (cv_equals_curr) {
      // m_cv overlaps with the curve of he_curr:
      m_intersect_he = he_curr;
      return true;
    }
    else if (cv_equals_next) {
      // m_cv overlaps with the curve of he_next:
      m_intersect_he = he_next;
      return true;
    }

    if (is_between) {
      // We can conclude that m_cv should be placed between he_curr and
      // he_next (in a clockwise order), and no overlap occurs.
      he = he_curr;
      return false;
    }

    // Proceed to the next halfedge around the vertex.
    he_curr = he_next;
    ++he_next;
  } while (he_curr != he_first);

  CGAL_error();                 // we should never reach here:
  return false;
}

//-----------------------------------------------------------------------------
// Direct the halfedge for the location of the given subcurve around a split
// point that occurs in the interior of a given edge, when the subcurve lies
// to the right of the split point.
//
template <typename Arrangement, typename ZoneVisitor>
typename Arrangement_zone_2<Arrangement, ZoneVisitor>::Halfedge_handle
Arrangement_zone_2<Arrangement, ZoneVisitor>::
_direct_intersecting_edge_to_right(const X_monotone_curve_2& cv_ins,
                                   const Point_2& cv_left_pt,
                                   Halfedge_handle query_he)
{
  // Make sure that the left endpoint of cv_ins lies on query_he.
  CGAL_exactness_assertion(m_geom_traits->compare_y_at_x_2_object()
                           (cv_left_pt, query_he->curve()) == EQUAL);

  // Check whether the given halfedge is directed to the right.
  const bool query_he_directed_right =
    (query_he->direction() == ARR_LEFT_TO_RIGHT);

  // Check whether the curve lies above of below the edge immediately to
  // the right of its left endpoint.
  const Comparison_result  pos_res =
    m_geom_traits->compare_y_at_x_right_2_object()(cv_ins, query_he->curve(),
                                                   cv_left_pt);

  if (pos_res == SMALLER) {
    // If m_cv below the curve associated with query_he, the relevant halfedge
    // is the one directed from right to left.
    if (query_he_directed_right) return (query_he->twin());
    else return (query_he);
  }
  else if (pos_res == LARGER) {
    // If m_cv below the curve associated with hh, the relevant halfedge
    // is the one directed from left to right.
    if (query_he_directed_right) return (query_he);
    else return (query_he->twin());
  }

  // The two curves are equal to the right of the left endpoint, so we have
  // an overlap.
  m_found_overlap = true;
  m_intersect_he = query_he;

  return (query_he);
}

//-----------------------------------------------------------------------------
// Direct the halfedge for the location of the given subcurve around a split
// point that occurs in the interior of a given edge, when the subcurve lies
// to the left of the split point.
//
template <typename Arrangement, typename ZoneVisitor>
typename Arrangement_zone_2<Arrangement, ZoneVisitor>::Halfedge_handle
Arrangement_zone_2<Arrangement, ZoneVisitor>::
_direct_intersecting_edge_to_left(const X_monotone_curve_2& cv_ins,
                                  Halfedge_handle query_he)
{
  // Make sure that the right endpoint of cv_ins lies on query_he.
  CGAL_exactness_assertion
    (m_geom_traits->compare_y_at_x_2_object()
     (m_geom_traits->construct_max_vertex_2_object()(cv_ins),
      query_he->curve()) == EQUAL);

  // Check whether the given halfedge is directed to the right.
  const bool query_he_directed_right =
    (query_he->direction() == ARR_LEFT_TO_RIGHT);

  // Check whether the curve lies above of below the edge (we use the curve
  // position predicate, as we know they cruves do not overlap and intersect
  // only at the split point).
  Comparison_result        pos_res =
      m_geom_traits->compare_y_position_2_object()(cv_ins, query_he->curve());

  if (pos_res == EQUAL) {
    // This can happen only when both endpoints of cv_ins lie on query_he,
    // for example (the ^-shaped polyline is associated with query_he and
    // the horizontal segment is cv_ins):
    //
    //      /\        .
    //     /  \       .
    //    +----+      .
    //   /      \     .
    //
    // In this case, we got a wrong result from compare_y_position(), as we
    // abused this predicate (since the two curves are not supposed to
    // intersect), so we now simply have to compare the two curves to the right
    // of cv_ins' left endpoint.
    pos_res = m_geom_traits->compare_y_at_x_right_2_object()
      (cv_ins, query_he->curve(),
       m_geom_traits->construct_min_vertex_2_object()(cv_ins));
  }

  if (pos_res == SMALLER) {
    // If cv_ins lies below the curve associated with query_he, we should
    // take the halfedge directed from right to left, so if query_he is
    // directed to the right, we return it twin.
    if (query_he_directed_right) return (query_he->twin());
    else return (query_he);
  }
  else {
    CGAL_assertion(pos_res != EQUAL);

    // If cv_ins lies above the curve associated with query_he, we should
    // take the halfedge directed from left to right, so if query_he is
    // directed to the left, we return it twin.
    if (! query_he_directed_right) return (query_he->twin());
    else return (query_he);
  }
}

//-----------------------------------------------------------------------------
// Get the next intersection of cv with the given halfedge.
//
template <typename Arrangement, typename ZoneVisitor>
typename Arrangement_zone_2<Arrangement, ZoneVisitor>::Optional_intersection
Arrangement_zone_2<Arrangement, ZoneVisitor>::
_compute_next_intersection(Halfedge_handle he,
                           bool skip_first_point,
                           bool& intersection_on_right_boundary)
{
  // Get a pointer to the curve associated with the halfedge.
  const X_monotone_curve_2* p_curve = &(he->curve());

  // Try to locate the intersections with this curve in the intersections map.
  Intersect_map_iterator iter = m_inter_map.find(p_curve);
  const Intersection_point* ip;
  const X_monotone_curve_2* icv;
  bool valid_intersection;

  intersection_on_right_boundary = false;
  if (iter != m_inter_map.end()) {
    // The intersections with the curve have already been computed.
    // Retrieve the intersections list from the map.
    Intersect_list& inter_list = iter->second;

    if (inter_list.empty()) return Optional_intersection();

    // Locate the first intersection that lies to the right of m_left_pt
    // (if the left point exists).
    while (! inter_list.empty()) {
      // Compare that current object with m_left_pt (if exists).
      ip = boost::get<Intersection_point>(&(inter_list.front()));

      if (m_left_on_boundary) {
        // The left end lie on the left boundary, so all intersections are
        // valid, as they lie to its right.
        valid_intersection = true;
      }
      else if (ip != nullptr) {
        if (m_has_right_pt && m_right_on_boundary &&
            m_geom_traits->equal_2_object()(ip->first, m_right_pt))
        {
          valid_intersection = true;
          intersection_on_right_boundary = true;
        }
        else {
          // We have a simple intersection point - make sure it lies to the
          // right of m_left_pt.
          valid_intersection =
            (m_geom_traits->compare_xy_2_object()(ip->first, m_left_pt) ==
             LARGER);
        }
      }
      else {
        // We have an overlapping subcurve.
        icv = boost::get<X_monotone_curve_2>(&(inter_list.front()));
        CGAL_assertion(icv != nullptr);

        if (m_geom_traits->is_closed_2_object()(*icv, ARR_MIN_END)) {
          // The curve has a valid left point - make sure it lies to the
          // right of m_left_pt.
          valid_intersection =
            (m_geom_traits->compare_xy_2_object()
             (m_geom_traits->construct_min_vertex_2_object()(*icv), m_left_pt) !=
             SMALLER);
        }
        else {
          // In this case the overlap is not valid.
          valid_intersection = false;
        }
      }

      // Found an intersection to m_left_pt's right.
      if (valid_intersection) return Optional_intersection(inter_list.front());

      // Discard the current intersection, which lies to m_left_pt's left.
      inter_list.pop_front();
    }

    // If we reached here, the list of intersections is empty:
    return Optional_intersection();
  }

  // The intersections with the curve have not been computed yet, so we
  // have to compute them now. Note that the first curve we intersect is
  // always the subcurve associated with the given halfegde and the second
  // curve is the one we insert. Even though the order seems unimportant, we
  // exploit this fact in some of the traits classes in order to optimize
  // computations.
  Intersect_list inter_list;
  bool is_first = true;

  m_geom_traits->intersect_2_object()(he->curve(), m_cv,
                                      std::back_inserter(inter_list));

  // Discard all intersection lying to the left of m_left_pt (if exists).
  while (! inter_list.empty()) {
    // Compare that current object with m_left_pt (if exists).
    ip = boost::get<Intersection_point>(&(inter_list.front()));

    if (ip != nullptr) {
      // We have a simple intersection point - if we don't have to skip it,
      // make sure it lies to the right of m_left_pt (if m_left_pt is on the
      // left boundary, all points lie to it right).
      if (is_first && skip_first_point) valid_intersection = false;
      else if (m_left_on_boundary) valid_intersection = true;
      else if (m_has_right_pt && m_right_on_boundary &&
               m_geom_traits->equal_2_object()(ip->first, m_right_pt))
      {
        valid_intersection = true;
        intersection_on_right_boundary = true;
      }
      else {
        valid_intersection =
          (m_geom_traits->compare_xy_2_object()(ip->first, m_left_pt) == LARGER);
      }
    }
    else if (m_left_on_boundary) {
      // The left end is on the boundary, so all overlapping curves are valid,
      // as they lie to its right.
      valid_intersection = true;
    }
    else {
      // We have an overlapping subcurve.
      icv = boost::get<X_monotone_curve_2>(&(inter_list.front()));
      CGAL_assertion(icv != nullptr);

      if (m_geom_traits->is_closed_2_object()(*icv, ARR_MIN_END)) {
        // The curve has a valid left point - make sure it lies to the
        // right of m_left_pt.
        valid_intersection =
          (m_geom_traits->compare_xy_2_object()
           (m_geom_traits->construct_min_vertex_2_object()(*icv), m_left_pt) !=
           SMALLER);
      }
      // In this case the overlap is not valid.
      else valid_intersection = false;
    }
    is_first = false;

    // Break, if an intersection to m_left_pt's right has been found.
    if (valid_intersection) break;

    // Discard the current intersection, which lies to m_left_pt's left.
    inter_list.pop_front();
  }

  // Insert the list of valid intersections into the map.
  m_inter_map[p_curve] = inter_list;

  // Return the first intersection object computed (may be empty).
  if (inter_list.empty()) return Optional_intersection();
  else return Optional_intersection(inter_list.front());
}

//-----------------------------------------------------------------------------
// Remove the next intersection of m_cv with the given halfedge from the map.
//
template <typename Arrangement, typename ZoneVisitor>
void Arrangement_zone_2<Arrangement, ZoneVisitor>::
_remove_next_intersection(Halfedge_handle he)
{
  // Get a pointer to the curve associated with the halfedge.
  const X_monotone_curve_2* p_curve = &(he->curve());

  // Locate the intersections with this curve in the intersections map.
  Intersect_map_iterator iter = m_inter_map.find(p_curve);

  CGAL_assertion(iter != m_inter_map.end());
  CGAL_assertion(! iter->second.empty());

  // Remove the first object in the list of intersections.
  iter->second.pop_front();
}

//-----------------------------------------------------------------------------
// Check if the given point lies completely to the left of the given egde.
//
template <typename Arrangement, typename ZoneVisitor>
bool Arrangement_zone_2<Arrangement, ZoneVisitor>::
_is_to_left_impl(const Point_2& p, Halfedge_handle he,
                 Arr_not_all_sides_oblivious_tag) const
{
  // Check the boundary conditions of the minimal end of the curve associated
  // with the given halfedge.
  const Arr_parameter_space ps_x =
    m_geom_traits->parameter_space_in_x_2_object()(he->curve(), ARR_MIN_END);

  // The minimal end of the curve is to the left of any other point:
  if (ps_x == ARR_LEFT_BOUNDARY) return false;

  const Arr_parameter_space ps_y =
    m_geom_traits->parameter_space_in_y_2_object()(he->curve(), ARR_MIN_END);

  if (ps_y != ARR_INTERIOR) {
    // Check if p is to the left of the minimal curve-end:
    const Comparison_result res =
      m_geom_traits->compare_x_point_curve_end_2_object()(p, he->curve(),
                                                          ARR_MIN_END);

    return ((res == SMALLER) || (res == EQUAL && ps_y == ARR_TOP_BOUNDARY));
  }

  // In case the minimal curve-end does not have boundary conditions, simply
  // compare p with the left endpoint of the curve.
  Vertex_const_handle v_left =
    (he->direction() == ARR_LEFT_TO_RIGHT) ? he->source() : he->target();

  return (m_geom_traits->compare_xy_2_object()(p, v_left->point()) == SMALLER);
}

//-----------------------------------------------------------------------------
// Check if the given point lies completely to the right of the given egde.
//
template <typename Arrangement, typename ZoneVisitor>
bool Arrangement_zone_2<Arrangement, ZoneVisitor>::
_is_to_right_impl(const Point_2& p, Halfedge_handle he,
                  Arr_not_all_sides_oblivious_tag) const
{
  // Check the boundary conditions of the maximal end of the curve associated
  // with the given halfedge.
  const Arr_parameter_space ps_x =
    m_geom_traits->parameter_space_in_x_2_object()(he->curve(), ARR_MAX_END);

    // The maximal end of the curve is to the right of any other point:
  if (ps_x == ARR_RIGHT_BOUNDARY) return false;

  const Arr_parameter_space ps_y =
    m_geom_traits->parameter_space_in_y_2_object()(he->curve(), ARR_MAX_END);

  if (ps_y != ARR_INTERIOR) {
    // Check if p is to the right of the maximal curve-end:
    const Comparison_result   res =
      m_geom_traits->compare_x_point_curve_end_2_object()(p, he->curve(),
                                                          ARR_MAX_END);

    return ((res == LARGER) || (res == EQUAL && ps_y == ARR_BOTTOM_BOUNDARY));
  }

  // In case the maximal curve-end does not have boundary conditions, simply
  // compare p with the right endpoint of the curve.
  Vertex_const_handle   v_right =
    (he->direction() == ARR_LEFT_TO_RIGHT) ? he->target() : he->source();

  return (m_geom_traits->compare_xy_2_object()(p, v_right->point()) == LARGER);
}

//-----------------------------------------------------------------------------
// Compute the (lexicographically) leftmost intersection of the query
// curve with a given halfedge on the boundary of a face in the arrangement.
template <typename Arrangement, typename ZoneVisitor>
void Arrangement_zone_2<Arrangement, ZoneVisitor>::
_leftmost_intersection(Ccb_halfedge_circulator he_curr, bool on_boundary,
                       bool& leftmost_on_right_boundary)
{
  // Obtain some geometry-traits functors.
  typename Traits_adaptor_2::Compare_xy_2 compare_xy =
    m_geom_traits->compare_xy_2_object();
  typename Traits_adaptor_2::Is_in_x_range_2 is_in_x_range =
    m_geom_traits->is_in_x_range_2_object();
  typename Traits_adaptor_2::Construct_min_vertex_2 min_vertex =
    m_geom_traits->construct_min_vertex_2_object();

  // If this edge is fictitious, skip it.
  if (he_curr->is_fictitious()) return;

  // If we have already found an intersection with the twin halfedge, we do not
  // have to compute intersections with the current halfedge.
  if (m_found_intersect && (m_intersect_he == he_curr->twin())) return;

  // If we already have an intersection point, compare it to the endpoints of
  // the curve associated with the current halfedge, in order to filter
  // unnecessary intersection computations. If the current x-monotone curve lies
  // entirely to the right of m_intersect_p, its intersection with m_cv (if any)
  // cannot lie to the left of this point. We therefore do not need to compute
  // this intersection.
  if (m_found_intersect && ! leftmost_on_right_boundary &&
      _is_to_left(m_intersect_p, he_curr))
    return;

  bool left_equals_curr_endpoint = false;
  if (on_boundary) {
    // Check if the left endpoint of the inserted curve (which is located on the
    // boundary of our face) equals one of the endpoints of the current
    // halfedge. If it equals the right endpoint of the current halfedge, we can
    // skip this edge, as there is no true overlap in the x-range. Otherwise, we
    // keep track of the fact that m_left_v is the left end-vertex of the
    // current halfedge.
    if (he_curr->target() == m_left_v) {
      left_equals_curr_endpoint = true;
      if (he_curr->direction() == ARR_LEFT_TO_RIGHT) return;
    }
    else if (he_curr->source() == m_left_v) {
      left_equals_curr_endpoint = true;
      if (he_curr->direction() == ARR_RIGHT_TO_LEFT) return;
    }
  }

  // Check whether the two curves overlap in their x-range (in order
  // to avoid unnecessary intersection computations).
  // In case there is no overlap, the two x-monotone curves obviously
  // do not intersect.
  if (! left_equals_curr_endpoint &&
      ((! m_left_on_boundary && _is_to_right(m_left_pt, he_curr)) ||
       ! is_in_x_range(m_cv, he_curr->curve())))
    return;

  // Compute the next intersection of m_cv and the current halfedge.
  bool intersection_on_right_boundary;
  Optional_intersection iobj =
    _compute_next_intersection(he_curr, left_equals_curr_endpoint,
                               intersection_on_right_boundary);

  if (iobj) {
    // We have found an intersection (either a simple point or an
    // overlapping x-monotone curve).
    const Intersection_point* int_p = boost::get<Intersection_point>(&*iobj);
    if (int_p != nullptr) {
      Point_2 ip = int_p->first;

      // Found a simple intersection point. Check if it is the leftmost
      // intersection point so far.
      if (! m_found_intersect ||
          (! intersection_on_right_boundary &&
           (leftmost_on_right_boundary ||
            compare_xy(ip, m_intersect_p) == SMALLER)))
      {
        // Store the leftmost intersection point and the halfedge handle.
        m_intersect_p = ip;
        m_ip_multiplicity = int_p->second;
        m_intersect_he = he_curr;
        m_found_overlap = false;
        leftmost_on_right_boundary = intersection_on_right_boundary;
      }
    }
    else {
      // We have located an overlapping curve. Assign ip as its left
      // endpoint.
      const X_monotone_curve_2* icv = boost::get<X_monotone_curve_2>(&*iobj);
      CGAL_assertion(icv != nullptr);
      Point_2 ip = min_vertex(*icv);

      // Check if this endpoint it is the leftmost intersection point so far.
      if (! m_found_intersect || compare_xy(ip, m_intersect_p) == SMALLER) {
        // Store the leftmost intersection point and the halfedge handle.
        m_intersect_p = ip;
        m_ip_multiplicity = 0;
        m_overlap_cv = *icv;
        m_intersect_he = he_curr;
        m_found_overlap = true;
      }
    }

    // Mark that we found an intersection.
    m_found_intersect = true;
  }
}

//-----------------------------------------------------------------------------
// Compute the (lexicographically) leftmost intersection of the query
// curve with the boundary of a given face in the arrangement.
//
template <typename Arrangement, typename ZoneVisitor>
void Arrangement_zone_2<Arrangement, ZoneVisitor>::
_leftmost_intersection_with_face_boundary(Face_handle face, bool on_boundary)
{
  // Mark that we have not found any intersection (or overlap) yet.
  m_found_intersect = false;
  m_found_overlap = false;
  m_found_iso_vert = false;

  // Obtain some geometry-traits functors.
  typename Traits_adaptor_2::Compare_xy_2 compare_xy =
    m_geom_traits->compare_xy_2_object();
  typename Traits_adaptor_2::Is_in_x_range_2 is_in_x_range =
    m_geom_traits->is_in_x_range_2_object();

  bool leftmost_on_right_boundary = false;

  // Traverse the face outer-boundaries; iterate through all outer CCBs.
  typedef typename Arrangement_2::Outer_ccb_iterator    Outer_ccb_iterator;
  for (Outer_ccb_iterator occb_it = face->outer_ccbs_begin();
       occb_it != face->outer_ccbs_end(); ++occb_it)
  {
    Ccb_halfedge_circulator he_first = *occb_it;
    Ccb_halfedge_circulator he_curr = he_first;
    do _leftmost_intersection(he_curr, on_boundary, leftmost_on_right_boundary);
    while (++he_curr != he_first);
  }

  // Traverse the face inner-boundaries; iterate through all inner CCBs (holes).
  typedef typename Arrangement_2::Inner_ccb_iterator    Inner_ccb_iterator;
  for (Inner_ccb_iterator iccb_it = face->inner_ccbs_begin();
       iccb_it != face->inner_ccbs_end(); ++iccb_it)
  {
    Ccb_halfedge_circulator he_first = *iccb_it;
    Ccb_halfedge_circulator he_curr = he_first;
    do _leftmost_intersection(he_curr, on_boundary, leftmost_on_right_boundary);
    while (++he_curr != he_first);
  }

  typename Traits_adaptor_2::Compare_y_at_x_2 compare_y_at_x =
    m_geom_traits->compare_y_at_x_2_object();

  // Traverse the isolated vertices inside the face (if there exist any), and
  // check whether an isolated vertex lies on the curve.
  typedef typename Arrangement_2::Isolated_vertex_iterator
    Isolated_vertex_iterator;
  for (Isolated_vertex_iterator iv_it = face->isolated_vertices_begin();
       iv_it != face->isolated_vertices_end(); ++iv_it)
  {
    // If the isolated vertex is not in the x-range of our curve, disregard it.
    if (! is_in_x_range(m_cv, iv_it->point())) continue;

    // If we already have an intersection point, compare it to the current
    // isolated vertex, in order to filter unnecessary computations.
    if (m_found_intersect && compare_xy(iv_it->point(), m_intersect_p) == LARGER)
      continue;

    // In case the isolated vertex lies on the curve, update the intersection
    // point accordingly.
    if ((compare_y_at_x(iv_it->point(), m_cv) == EQUAL) &&
        (! m_has_left_pt || (compare_xy(iv_it->point(), m_left_pt) == LARGER)))
    {
      m_intersect_v = iv_it;
      m_intersect_p = m_intersect_v->point();
      m_ip_multiplicity = 0;
      m_found_intersect = true;
      m_found_iso_vert = true;
    }
  } // End:: traversal of the isolated vertices inside the face.

  // Remove the next intersection associated with m_intersect_he, as we have
  // now reported it and do not want to encounter it again.
  if (m_found_intersect && !m_found_iso_vert)
    _remove_next_intersection(m_intersect_he);
}

//-----------------------------------------------------------------------------
// Compute the zone of an x-monotone curve in a given arrangement face.
// The left endpoint of the curve either lies in the face interior or on
// the boundary of the face.
//
template <typename Arrangement, typename ZoneVisitor>
bool Arrangement_zone_2<Arrangement, ZoneVisitor>::
_zone_in_face(Face_handle face, bool on_boundary)
{
  // Obtain some geometry-traits functors.
  typename Traits_adaptor_2::Equal_2 equal = m_geom_traits->equal_2_object();

  CGAL_precondition((! on_boundary &&
                     (((m_left_v == m_invalid_v) &&
                       (m_left_he == m_invalid_he)) ||
                      m_left_v->is_isolated())) ||
                    (on_boundary && (m_left_he != m_invalid_he)));

  // Find the first intersection of the curve with the face boundary.
  _leftmost_intersection_with_face_boundary(face, on_boundary);

  if (! m_found_intersect) {
    // Notify the visitor that the entire curve lies within the given face,
    // such that its right endpoint is not incident to any arrangement feature.
    m_visitor->found_subcurve(m_cv, face, m_left_v, m_left_he,
                              m_invalid_v, m_invalid_he);

    // Inidicate that we are done with the zone-computation process.
    return true;
  }

  // In this case m_found_intersect is true and m_intersect_he is the edge that
  // m_cv next intersects (or overlaps). If m_found_overlap is also true,
  // then m_overlap_cv is set and m_intersect_p is the left endpoint of the
  // overlapping subcurve. Otherwise, m_intersect_p is a simple intersection
  // point.
  // Alternatively, if m_found_iso_vert is true, then the next intersection point
  // m_intersect_p lies on the isolated vertex m_intersect_v.
  bool done = false;

  if (m_has_right_pt && equal(m_intersect_p, m_right_pt)) {
    // If the intersection point is m_cv's right endpoint, the interior of cv
    // does not intersect any existing halfedge. In this case, we only have
    // to insert m_cv to the arrangement and we are done.
    m_sub_cv1 = m_cv;
    done = true;
  }
  else {
    // Split cv at the intersection point.
    m_geom_traits->split_2_object()(m_cv, m_intersect_p, m_sub_cv1, m_sub_cv2);

    // Set m_cv to be the remaining portion.
    m_has_left_pt = true;
    m_left_on_boundary = false;
    m_left_pt = m_intersect_p;
    m_cv = m_sub_cv2;
  }

  const X_monotone_curve_2*p_orig_curve = nullptr;

  if (! m_found_iso_vert) {
    // Check whether m_intersect_p coincides with one of the end-vertices of the
    // halfedge that m_cv intersects.
    if (! m_intersect_he->source()->is_at_open_boundary() &&
        equal(m_intersect_p, m_intersect_he->source()->point()))
    {
      // We know that the right endpoint of m_sub_cv1 lies on the source vertex:
      m_right_v = m_intersect_he->source();
      m_right_he = m_invalid_he;
    }
    else if (! m_intersect_he->target()->is_at_open_boundary() &&
             equal(m_intersect_p, m_intersect_he->target()->point()))
    {
      // We know that the right endpoint of m_sub_cv1 lies on the target vertex:
      m_right_v = m_intersect_he->target();
      m_right_he = m_invalid_he;
    }
    else {
      // The right endpoint of m_sub_cv1 lies on the interior of m_intersect_he:
      // Obtain the halfedge with the correct direction (which should be the
      // predecessor of m_sub_cv1 if we split the edge around this vertex).
      m_right_v = m_invalid_v;
      m_right_he = _direct_intersecting_edge_to_left(m_sub_cv1, m_intersect_he);
    }

    // Store the curve currently associated with the intersecting halfedge.
    p_orig_curve = &(m_intersect_he->curve());
  }
  else {
    // The right endpoint of the subcurve coincides with an isolated vertex:
    m_right_v = m_intersect_v;
    m_right_he = m_invalid_he;
  }
  // Notify the visitor that the left endpoint of the first subcurve is
  // located within the current face and both its endpoint are located
  // on its boundary.
  Visitor_result visitor_res = m_visitor->found_subcurve(m_sub_cv1, face,
                                                         m_left_v, m_left_he,
                                                         m_right_v, m_right_he);

  // Check if we are done (either we have no remaining curve or if the
  // visitor has indicated we should end the process).
  if (done || visitor_res.second) return true;

  // Move to the remaining portion of the curve, whose left endpoint is the
  // same as the right endpoint of m_sub_cv1. Note that we check if the visitor
  // has inserted the subcurve (in which case it should return a handle to
  // the resulting halfedge).
  Halfedge_handle inserted_he = visitor_res.first;

  if (inserted_he != m_invalid_he) {
    if (m_right_v == m_invalid_v) {
      // If the right endpoint of the subcurve we have just detected was
      // not associated with an existing vertex, the inserted halfedge is
      // now targeted toward a newly created vertex that splits m_intersect_he
      // into two halfedges: (a) the next halfedge after inserted_he and (b)
      // the previous halfedge before inserted_he's twin.
      // The two halfedges (a) and (b) are now associated with the two
      // subcurves that result from splitting m_intersect_he->curve() at the
      // intersection point we have just detected, one extends to the left
      // and one to the right of this split point.
      const X_monotone_curve_2* p_left_subcurve = nullptr;
      const X_monotone_curve_2* p_right_subcurve = nullptr;

      if (inserted_he->next()->direction() == ARR_LEFT_TO_RIGHT) {
        // The next halfedge extends to the right of the split point:
        p_left_subcurve = &(inserted_he->twin()->prev()->curve());
        p_right_subcurve = &(inserted_he->next()->curve());
      }
      else {
        // The next halfedge extends to the left of the split point:
        p_right_subcurve = &(inserted_he->twin()->prev()->curve());
        p_left_subcurve = &(inserted_he->next()->curve());
      }

      // Associate the intersection list of the original curve with the
      // right subcurve, while we can associate an empty list with the
      // left subcurve, as we are now done with it.
      Intersect_map_iterator iter = m_inter_map.find(p_orig_curve);
      Intersect_list empty_inter_list;

      m_inter_map[p_right_subcurve] = iter->second;
      m_inter_map[p_left_subcurve] = empty_inter_list;

      // If necessary, erase the original curve from the intersection map.
      if (p_orig_curve != p_right_subcurve && p_orig_curve !=p_left_subcurve)
        m_inter_map.erase(p_orig_curve);
    }

    if (m_found_overlap && m_right_v == m_invalid_v) {
      // In case we have split the overlapping m_intersect_he, it now refers
      // to the wrong halfedge. the overlapping edge is either the successor
      // of the inserted halfedge or the predecessor of its twin, depending
      // on which one of these halfedges lies to the right of the split point.
      if (inserted_he->next()->direction() == ARR_LEFT_TO_RIGHT)
        // The successor is directed to the right:
        m_intersect_he = inserted_he->next();
      else {
        // The predecessor is directed to the left:
        CGAL_assertion(inserted_he->twin()->prev()->direction() ==
                       ARR_RIGHT_TO_LEFT);

        m_intersect_he = inserted_he->twin()->prev();
      }
    }

    // The visitor has created an edge that corresponds to m_sub_cv1 and
    // inserted it into the arrangement. In this case, m_left_pt should be
    // associated with the target vertex of the new halfedge.
    CGAL_assertion(equal(m_left_pt, inserted_he->target()->point()));

    m_left_v = inserted_he->target();

    // If m_right_he is known, it is possible to set m_left_he according to the
    // geometric information we have.
    if (m_right_he != m_invalid_he) {
      if ((m_ip_multiplicity % 2) == 1) {
        // m_cv crosses m_right_he (which is now split into two), so the
        // remaining portion must be inserted after the next halfedge going
        // clockwise around m_left_v:
        //
        //              \   .                            .
        //               \ . remaining portion of m_cv   .
        //                x                              .
        //   inserted_he / \                             .
        //              /   \                            .
        m_left_he = inserted_he->next()->twin();
      }
      else if (m_ip_multiplicity != 0)
        // We have a tangency point. If m_right_he is directed from left to
        // right, we take the inserted halfedge to be m_left_he, otherwise
        // m_right_he itself becomes m_left_he:
        m_left_he = (m_right_he->direction() == ARR_LEFT_TO_RIGHT) ?
          inserted_he : m_right_he;
      else {
        // Mutliplicity is unkown:
        m_left_he = m_invalid_he;
      }
    }
    else {
      // In case m_left_v used to be an isolated vertex, we know that the
      // inserted halfedge is its only incident halfedge and we can use it.
      // Otherwise, we do not know the identity of m_left_he.
      m_left_he = (m_found_iso_vert) ? inserted_he : m_invalid_he;
    }
  }
  else {
    // The visitor has not created a new edge. We proceed using the previously
    // computed arrangement features.
    m_left_v = m_right_v;

    if (m_right_he != m_invalid_he) {
      // In case m_cv crosses the interior of the m_right_he halfedge (the
      // multiplicity of the intersection is odd), we know that the ramaining
      // portion of the curve lies in the face incident to the twin halfedge.
      // If the multiplicity is known and is even, we stay with the same
      // halfedge.
      if ((m_ip_multiplicity % 2) == 1) m_left_he = m_right_he->twin();
      else if (m_ip_multiplicity != 0) m_left_he = m_right_he;
      else m_left_he = m_invalid_he;
    }
    else m_left_he = m_invalid_he;
  }

  // We are not done with the zone-computation process yet:
  return false;
}

//-----------------------------------------------------------------------------
// Compute the zone of an overlapping subcurve m_overlap_cv of m_cv and the
// curve currently associated with m_intersect_he.
//
template <typename Arrangement, typename ZoneVisitor>
bool Arrangement_zone_2<Arrangement, ZoneVisitor>::_zone_in_overlap()
{
  // Obtain some geometry-traits functors.
  typename Traits_adaptor_2::Equal_2 equal = m_geom_traits->equal_2_object();

  // Check if the right end of m_overlap_cv is bounded. If so, compute its
  // right endpoint.
  const bool cv_has_right_pt =
    m_geom_traits->is_closed_2_object()(m_overlap_cv, ARR_MAX_END);

  Point_2 cv_right_pt;

  if (cv_has_right_pt)
    cv_right_pt = m_geom_traits->construct_max_vertex_2_object()(m_overlap_cv);

  // Get right end-vertex of the overlapping halfedge m_intersect_he. Also make
  // sure that the overlapping halfedge is always directed to the right.
  Vertex_handle he_right_v;

  if (m_intersect_he->direction() == ARR_LEFT_TO_RIGHT) {
    he_right_v = m_intersect_he->target();
  }
  else {
    he_right_v = m_intersect_he->source();
    m_intersect_he = m_intersect_he->twin();
  }

  // Compare the two right endpoints. Note that m_overlap_cv cannot extend to
  // the right longer than the halfedge it overlaps. Thus, if the curve is not
  // bounded, the right vertex of m_intersect_he must lie on open boundary as
  // well.
  if (! cv_has_right_pt) {
    CGAL_assertion_code
      (const Arr_parameter_space cv_ps_x =
       m_geom_traits->parameter_space_in_x_2_object()(m_overlap_cv, ARR_MAX_END);
       const Arr_parameter_space cv_ps_y =
       m_geom_traits->parameter_space_in_y_2_object()(m_overlap_cv, ARR_MAX_END);
       );
    CGAL_assertion(he_right_v->parameter_space_in_x() == cv_ps_x &&
                   he_right_v->parameter_space_in_y() == cv_ps_y);

    m_right_v = he_right_v;
  }
  else {
    // In this case m_overlap_cv has a finite right endpoint. In this case,
    // if the right vertex of m_intersect_he is associated with a finite point,
    // we check whether it is equal to cv_right_pt. Otherwise, we know that
    // m_intersect_he extends to the the right of m_overlap_cv, and there is no
    // vertex currently associated with m_overlap_cv's right endpoint.
    if (! he_right_v->is_at_open_boundary() &&
        equal(cv_right_pt, he_right_v->point()))
    {
      // The overlap is with the entire halfedge. In this case we set the
      // right end-vertex of the overlapping zone.
      m_right_v = he_right_v;
    }
    else {
      // In this case m_intersect_he overlaps just a portion of prev_he.
      // The right end-vertex of the overlapping zone is not known.
      m_right_v = m_invalid_v;
    }
  }

  // Store the curve currently associated with the overlapping halfedge.
  const X_monotone_curve_2* p_orig_curve = &(m_intersect_he->curve());

  // Notify the visitor on the overlapping zone.
  Visitor_result visitor_res =
    m_visitor->found_overlap(m_overlap_cv, m_intersect_he, m_left_v, m_right_v);

  // If the visitor has indicated we should halt the process, or it the right
  // endpoint of the overlapping curve is the right endpoint of m_cv then we are
  // done (or both extend to an open boundary).
  if (visitor_res.second ||
      (cv_has_right_pt && m_has_right_pt && equal(cv_right_pt, m_right_pt)) ||
      (! cv_has_right_pt && ! m_has_right_pt))
  {
    return true;
  }

  // Erase the original curve from the intersection map, so we will have to
  // recompute intersections with it in the future.
  m_inter_map.erase(p_orig_curve);

  // Mark that we have dealt with the overlap.
  m_found_overlap = false;

  // Split m_cv at right endpoint of the overlapping curve.
  m_geom_traits->split_2_object()(m_cv, cv_right_pt, m_sub_cv1, m_sub_cv2);

  // Set m_cv to be the remaining portion.
  m_has_left_pt = true;
  m_left_on_boundary = false;
  m_left_pt = cv_right_pt;
  m_cv = m_sub_cv2;

  // Move to the remaining portion of the curve, whose left endpoint is the
  // same as the right endpoint of the overlapping curve. Note that we check
  // if the visitor has inserted the subcurve (in which case it should return
  // a handle to the resulting halfedge).
  Halfedge_handle updated_he = visitor_res.first;

  if (updated_he != m_invalid_he) {
    // In this case, m_left_pt should be associated with the target vertex of
    // the updated halfedge.
    CGAL_assertion(equal(m_left_pt, updated_he->target()->point()));

    m_left_v = updated_he->target();
  }
  else m_left_v = m_right_v;

  m_left_he = m_invalid_he;

  // We are not done with the zone-computation process yet:
  return false;
}

} //namespace CGAL

#endif
