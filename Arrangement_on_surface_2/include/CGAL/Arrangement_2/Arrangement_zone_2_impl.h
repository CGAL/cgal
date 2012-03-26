// Copyright (c) 2005,2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Ron Wein          <wein@post.tau.ac.il>
//                 Efi Fogel         <efif@post.tau.ac.il>
//                 (based on old version by Eyal Flato)

#ifndef CGAL_ARRANGEMENT_ZONE_2_IMPL_H
#define CGAL_ARRANGEMENT_ZONE_2_IMPL_H

/*! \file
 * Member-function definitions for the Arrangement_zone_2 class.
 */

namespace CGAL {

//-----------------------------------------------------------------------------
// Initialize the zone-computation process with a given curve and an object
// that wraps the location of the curve's left end.
//
template<class Arrangement, class ZoneVisitor>
void Arrangement_zone_2<Arrangement,ZoneVisitor>::
init_with_hint(const X_monotone_curve_2& _cv, const Object& _obj)
{
  // Set the curve and check whether its ends are bounded, therefore
  // associated with valid endpoints.
  cv = _cv;
  
  if (m_geom_traits->is_closed_2_object()(cv, ARR_MIN_END))
  {
    // The left endpoint is valid.
    const Arr_parameter_space  ps_x1 =
      m_geom_traits->parameter_space_in_x_2_object()(cv, ARR_MIN_END);
    const Arr_parameter_space  ps_y1 =
      m_geom_traits->parameter_space_in_y_2_object()(cv, ARR_MIN_END);
    has_left_pt = true;
    left_on_boundary = (ps_x1 != ARR_INTERIOR || ps_y1 != ARR_INTERIOR);
    left_pt = m_geom_traits->construct_min_vertex_2_object() (cv);
  }
  else
  {
    // The left end of the curve lies on open boundary.
    has_left_pt = false;
    left_on_boundary = true;
  }
  
  if (m_geom_traits->is_closed_2_object()(cv, ARR_MAX_END))
  {
    // The right endpoint is valid.
    const Arr_parameter_space  ps_x2 =
      m_geom_traits->parameter_space_in_x_2_object()(cv, ARR_MAX_END);
    const Arr_parameter_space  ps_y2 =
      m_geom_traits->parameter_space_in_y_2_object()(cv, ARR_MAX_END);
    has_right_pt = true;
    right_on_boundary = (ps_x2 != ARR_INTERIOR || ps_y2 != ARR_INTERIOR);
    right_pt = m_geom_traits->construct_max_vertex_2_object() (cv);
  }
  else
  {
    // The right end of the curve lies on open boundary.
    has_right_pt = false;
    right_on_boundary = true;
  }
  
  // Set the object that represents the location of the left end of the curve
    // in the arrangement.
  obj = _obj;

  return;
}

//-----------------------------------------------------------------------------
// Compute the zone of the given curve and issue the apporpriate
// notifications for the visitor.
//
template<class Arrangement, class ZoneVisitor>
void Arrangement_zone_2<Arrangement,ZoneVisitor>::compute_zone ()
{
  // Initialize flags and set all handles to be invalid.
  bool    done = false;

  found_intersect = false;
  found_overlap = false;
  found_iso_vert = false;

  left_v = invalid_v;
  left_he = invalid_he;
  right_v = invalid_v;
  right_he = invalid_he;

  // Locate the arrangement feature containing the left endpoint of the
  // curve (currently obj stores the object containing it).
  const Vertex_const_handle    *vh;
  const Halfedge_const_handle  *hh;
  const Face_const_handle      *fh;

  if ((vh = object_cast<Vertex_const_handle>(&obj)) != NULL)
  {
    CGAL_assertion (has_left_pt);

    // The left endpoint coincides with an existing vertex:
    left_v = arr.non_const_handle (*vh);

    if (left_on_boundary)
    {
      // Use the accessor to locate the predecessor edge, in case the left
      // endpoint has boundary conditions.
      const Arr_parameter_space  ps_x =
        m_geom_traits->parameter_space_in_x_2_object()(cv, ARR_MIN_END);
      const Arr_parameter_space  ps_y =
        m_geom_traits->parameter_space_in_y_2_object()(cv, ARR_MIN_END);

      left_he = arr_access.locate_around_boundary_vertex (left_v,
                                                          cv, ARR_MIN_END,
                                                          ps_x, ps_y);
    }
  }
  else if ((hh = object_cast<Halfedge_const_handle>(&obj)) != NULL)
  {
    if (has_left_pt)
    {
      // Obtain the right halfedge from the halfedge-pair containing left_pt
      // in their interior.
      left_he = _direct_intersecting_edge_to_right (cv, left_pt,
                                                    arr.non_const_handle(*hh));

      // Handle overlaps.
      if (found_overlap)
      {
        // In this case cv overlaps the curve associated with intersect_he.
        // Compute the overlapping subcurve.
        bool    dummy;
        obj = _compute_next_intersection (intersect_he, false, dummy);
        
        overlap_cv = object_cast<X_monotone_curve_2> (obj);
        
        // Remove the overlap from the map.
        _remove_next_intersection (intersect_he);
        
        // Compute the overlap zone.
        done = _zone_in_overlap ();
      }
    }
    else
    {
      // In case the unbounded left end conicides with an edge, then our curve
      // overlaps the curve associated with this edge.
      intersect_he = arr.non_const_handle (*hh);

      bool  dummy;
      obj = _compute_next_intersection (intersect_he, false, dummy);
        
      overlap_cv = object_cast<X_monotone_curve_2> (obj);
        
      // Remove the overlap from the map.
      _remove_next_intersection (intersect_he);
        
      // Compute the overlap zone.
      done = _zone_in_overlap ();
    }
  }
  else
  {
    // The left endpoint lies inside a face.
    fh = object_cast<Face_const_handle>(&obj);

    CGAL_assertion_msg(fh != NULL,
		       "Invalid object returned by the point-location query.");

    // Compute the zone of the curve at the interior of the face.
    done = _zone_in_face (arr.non_const_handle(*fh),
                          false);      // left_pt is not on the face boundary.

    // In case we have just discovered an overlap, compute the overlapping
    // zone as well.
    if (! done && found_overlap)
    {
      done = _zone_in_overlap ();
    }
  }

  // Compute the zone of the curve (or what is remaining of it) in the
  // arrangement, starting from the current position we have computed.
  while (! done)
  {
    // Check if we know the face the curve is going to penetrate now.
    if (left_he == invalid_he)
    {
      if (left_v != invalid_v)
      {
        // We know the vertex that coincides with the left endpoint of cv.
        if (! left_v->is_isolated())
        {
          // Locate the curve around the left_v vertex - that is, find a
          // halfedge left_he such that cv should be placed between left_he
          // and its current successor around the vertex, going in a clockwise
          // order.
          found_overlap = _find_prev_around_vertex (left_v, left_he);
        }
        else
        {
          // left_v is an isolated vertex.
          found_iso_vert = true;
        }
      }
      else
      {
        CGAL_assertion (right_he != invalid_he);

        // In this case right_he is the halfedge that the left portion of cv
        // intersected, and we obtain left_he by comparing the remaining
        // portion of cv with the curve associated with this edge.
        left_he = _direct_intersecting_edge_to_right (cv, left_pt, right_he);
      }

      if (found_overlap)
      {
        // In this case cv overlaps the curve associated with intersect_he.
        // Compute the overlapping subcurve to the right of curr_v.
        bool  dummy;
        obj = _compute_next_intersection (intersect_he, false, dummy);

        overlap_cv = object_cast<X_monotone_curve_2> (obj);

        // Remove the overlap from the map.
        _remove_next_intersection (intersect_he);

        // Compute the overlap zone and continue to the end of the loop.
        done = _zone_in_overlap ();
        continue;
      }
    }

    if (left_v == invalid_v || ! left_v->is_isolated())
    {
      // At this point we can compute the zone of cv starting from the left_he
      // inside its incident face.
      done = _zone_in_face (left_he->face(), true);
      // left_pt is on the face boundary.
    }
    else
    {
      // Compute the zone of cv starting from the face that contains the
      // isolated vertex left_v.
      done = _zone_in_face (left_v->face(), false);
      // left_pt is not on the face boundary.
    }

    // In case we have just discovered an overlap, compute the overlapping
    // zone as well.
    if (! done && found_overlap)
    {
      done = _zone_in_overlap();
    }
  }

  // Clear the intersections map.
  inter_map.clear();

  return;
}

//-----------------------------------------------------------------------------
// Find a face containing the query curve cv around the given vertex.
// In case an overlap occurs, sets intersect_he to be the overlapping edge.
//
template<class Arrangement, class ZoneVisitor>
bool Arrangement_zone_2<Arrangement,ZoneVisitor>::_find_prev_around_vertex
    (Vertex_handle v,
     Halfedge_handle& he)
{
  // Go over the incident halfedges of v, going in a clockwise order.
  typename Arrangement_2::Halfedge_around_vertex_circulator he_first;
  typename Arrangement_2::Halfedge_around_vertex_circulator he_curr;
  bool                                                      cv_equals_curr;
  typename Arrangement_2::Halfedge_around_vertex_circulator he_next;
  bool                                                      cv_equals_next;
  bool                                                      is_between;

  he_first = v->incident_halfedges();
  he_curr = he_first;
  he_next = he_curr;
  ++he_next;

  if (he_curr == he_next)
  {
    // In case there is just a single incident halfedge around v,
    // we should insert cv right after this halfedge.
    he = he_curr;

    // Note that cv extends to the right of v. In case the single
    // halfedge also extends to the right of v (its source is to
    // the right), check if an overlap occurs.
    if (he->direction() == ARR_RIGHT_TO_LEFT &&
        (m_geom_traits->compare_y_at_x_right_2_object() (he->curve(), cv,
                                                         v->point()) == EQUAL))
    {
      // Mark that an overlap occurs:
      intersect_he = he_curr;
      return (true);
    }

    // We have no overlap:
    return (false);
  }

  // Find the face containing cv around the vertex.
  typename Traits_adaptor_2::Is_between_cw_2  is_between_cw =
    m_geom_traits->is_between_cw_2_object();

  do
  {
    // Check if it is possible to insert cv in between the current curve
    // and the next curve, going in a clockwise direction around v.
    is_between = is_between_cw (cv, true,
                                he_curr->curve(), 
                                (he_curr->direction() == ARR_RIGHT_TO_LEFT),  
                                he_next->curve(),
                                (he_next->direction() == ARR_RIGHT_TO_LEFT),
                                v->point(),
                                cv_equals_curr, cv_equals_next);

    // Check the case of overlaps:
    if (cv_equals_curr)
    {
      // cv overlaps with the curve of he_curr:
      intersect_he = he_curr;
      return (true);
    }
    else if (cv_equals_next)
    {
      // cv overlaps with the curve of he_next:
      intersect_he = he_next;
      return (true);
    }

    if (is_between)
    {
      // We can conclude that cv should be placed between he_curr and
      // he_next (in a clockwise order), and no overlap occurs.
      he = he_curr;
      return (false);
    }

    // Proceed to the next halfedges around the vertex.
    ++he_curr;
    ++he_next;

  } while (he_curr != he_first);

  // We should never reach here:
  CGAL_error();
  return (false);
}

//-----------------------------------------------------------------------------
// Direct the halfedge for the location of the given subcurve around a split
// point that occurs in the interior of a given edge, when the subcurve lies
// to the right of the split point.
//
template<class Arrangement, class ZoneVisitor>
typename Arrangement_zone_2<Arrangement,ZoneVisitor>::Halfedge_handle
Arrangement_zone_2<Arrangement,ZoneVisitor>::_direct_intersecting_edge_to_right
    (const X_monotone_curve_2& cv_ins,
     const Point_2& cv_left_pt,
     Halfedge_handle query_he)
{
  // Make sure that the left endpoint of cv_ins lies on query_he.
  CGAL_exactness_assertion (m_geom_traits->compare_y_at_x_2_object()
                            (cv_left_pt, query_he->curve()) == EQUAL);

  // Check whether the given halfedge is directed to the right.
  const bool               query_he_directed_right = 
                                   (query_he->direction() == ARR_LEFT_TO_RIGHT);

  // Check whether the curve lies above of below the edge immediately to
  // the right of its left endpoint.
  const Comparison_result  pos_res =
    m_geom_traits->compare_y_at_x_right_2_object() (cv_ins, query_he->curve(),
                                                    cv_left_pt);

  if (pos_res == SMALLER)
  {
    // If cv below the curve associated with query_he, the relevant halfedge
    // is the one directed from right to left.
    if (query_he_directed_right)
      return (query_he->twin());
    else
      return (query_he);
  }
  else if (pos_res == LARGER)
  {
    // If cv below the curve associated with hh, the relevant halfedge
    // is the one directed from left to right.
    if (query_he_directed_right)
      return (query_he);
    else
      return (query_he->twin());
  }

  // The two curves are equal to the right of the left endpoint, so we have
  // an overlap.
  found_overlap = true;
  intersect_he = query_he;

  return (query_he);
}

//-----------------------------------------------------------------------------
// Direct the halfedge for the location of the given subcurve around a split
// point that occurs in the interior of a given edge, when the subcurve lies
// to the left of the split point.
//
template<class Arrangement, class ZoneVisitor>
typename Arrangement_zone_2<Arrangement,ZoneVisitor>::Halfedge_handle
Arrangement_zone_2<Arrangement,ZoneVisitor>::
_direct_intersecting_edge_to_left (const X_monotone_curve_2& cv_ins,
                                   Halfedge_handle query_he)
{
  // Make sure that the right endpoint of cv_ins lies on query_he.
  CGAL_exactness_assertion
    (m_geom_traits->compare_y_at_x_2_object()
     (m_geom_traits->construct_max_vertex_2_object()(cv_ins),
      query_he->curve()) == EQUAL);

  // Check whether the given halfedge is directed to the right.
  const bool               query_he_directed_right =
    (query_he->direction() == ARR_LEFT_TO_RIGHT);

  // Check whether the curve lies above of below the edge (we use the curve
  // position predicate, as we know they cruves do not overlap and intersect
  // only at the split point).
  Comparison_result        pos_res =
      m_geom_traits->compare_y_position_2_object() (cv_ins, query_he->curve());

  if (pos_res == EQUAL)
  {
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
       m_geom_traits->construct_min_vertex_2_object() (cv_ins));
  }

  if (pos_res == SMALLER)
  {
    // If cv_ins lies below the curve associated with query_he, we should
    // take the halfedge directed from right to left, so if query_he is
    // directed to the right, we return it twin.
    if (query_he_directed_right)
      return (query_he->twin());
    else
      return (query_he);
  }
  else
  {
    CGAL_assertion (pos_res != EQUAL);

    // If cv_ins lies above the curve associated with query_he, we should
    // take the halfedge directed from left to right, so if query_he is
    // directed to the left, we return it twin.
    if (! query_he_directed_right)
      return (query_he->twin());
    else
      return (query_he);
  }
}

//-----------------------------------------------------------------------------
// Get the next intersection of cv with the given halfedge.
//
template<class Arrangement, class ZoneVisitor>
CGAL::Object
Arrangement_zone_2<Arrangement,ZoneVisitor>::
_compute_next_intersection (Halfedge_handle he,
                            bool skip_first_point,
                            bool& intersection_on_right_boundary)
{
  // Get a pointer to the curve associated with the halfedge.
  const X_monotone_curve_2  *p_curve = &(he->curve());

  // Try to locate the intersections with this curve in the intersections map.
  Intersect_map_iterator    iter = inter_map.find (p_curve);
  const Intersect_point_2  *ip;
  const X_monotone_curve_2 *icv;
  bool                      valid_intersection;

  intersection_on_right_boundary = false;
  if (iter != inter_map.end())
  {
    // The intersections with the curve have already been computed.
    // Retrieve the intersections list from the map.
    Intersect_list&          inter_list = iter->second;

    if (inter_list.empty())
      return CGAL::Object();

    // Locate the first intersection that lies to the right of left_pt
    // (if the left point exists).
    while (! inter_list.empty())
    {
      // Compare that current object with left_pt (if exists).
      ip = object_cast<Intersect_point_2> (&(inter_list.front()));

      if (left_on_boundary)
      {
        // The left end lie on the left boundary, so all intersections are
        // valid, as they lie to its right.
        valid_intersection = true;
      }
      else if (ip != NULL)
      {
        if (has_right_pt && right_on_boundary &&
            m_geom_traits->equal_2_object() (ip->first, right_pt))
        {
          valid_intersection = true;
          intersection_on_right_boundary = true;
        }
        else
        {
          // We have a simple intersection point - make sure it lies to the
          // right of left_pt.
          valid_intersection =
            (m_geom_traits->compare_xy_2_object() (ip->first, left_pt) ==
             LARGER);
        }
      }
      else
      {
        // We have an overlapping subcurve.
        icv = object_cast<X_monotone_curve_2> (&(inter_list.front()));
        CGAL_assertion (icv != NULL);

        if (m_geom_traits->is_closed_2_object()(*icv, ARR_MIN_END))
        {
          // The curve has a valid left point - make sure it lies to the
          // right of left_pt.
          valid_intersection =
            (m_geom_traits->compare_xy_2_object()
             (m_geom_traits->construct_min_vertex_2_object()(*icv), left_pt) !=
             SMALLER);
        }
        else
        {
          // In this case the overlap is not valid.
          valid_intersection = false;
        }
      }

      if (valid_intersection)
        // Found an intersection to left_pt's right.
        return (inter_list.front());

      // Discard the current intersection, which lies to left_pt's left.
      inter_list.pop_front();
    }

    // If we reached here, the list of intersections is empty:
    return CGAL::Object();
  }

  // The intersections with the curve have not been computed yet, so we
  // have to compute them now. Note that the first curve we intersect is
  // always the subcurve associated with the given halfegde and the second
  // curve is the one we insert. Even though the order seems unimportant, we
  // exploit this fact in some of the traits classes in order to optimize
  // computations.
  Intersect_list           inter_list;
  bool                     is_first = true;

  m_geom_traits->intersect_2_object() (he->curve(), cv,
                                       std::back_inserter(inter_list));

  // Discard all intersection lying to the left of left_pt (if exists).
  while (! inter_list.empty())
  {
    // Compare that current object with left_pt (if exists).
    ip = object_cast<Intersect_point_2> (&(inter_list.front()));

    if (ip != NULL)
    {
      // We have a simple intersection point - if we don't have to skip it,
      // make sure it lies to the right of left_pt (if left_pt is on the left
      // boundary, all points lie to it right).
      if (is_first && skip_first_point)
      {
        valid_intersection = false;
      }
      else if (left_on_boundary)
      {
        valid_intersection = true;
      }
      else if (has_right_pt && right_on_boundary &&
               m_geom_traits->equal_2_object() (ip->first, right_pt))
      {
        valid_intersection = true;
        intersection_on_right_boundary = true;
      }
      else
      {
        valid_intersection =
          (m_geom_traits->compare_xy_2_object() (ip->first, left_pt) == LARGER);
      }
    }
    else if (left_on_boundary)
    {
      // The left end is on the boundary, so all overlapping curves are valid,
      // as they lie to its right.
      valid_intersection = true;
    }
    else
    {
      // We have an overlapping subcurve.
      icv = object_cast<X_monotone_curve_2> (&(inter_list.front()));
      CGAL_assertion (icv != NULL);

      if (m_geom_traits->is_closed_2_object() (*icv, ARR_MIN_END))
      {
        // The curve has a valid left point - make sure it lies to the
        // right of left_pt.
        valid_intersection =
          (m_geom_traits->compare_xy_2_object()
           (m_geom_traits->construct_min_vertex_2_object()(*icv),
            left_pt) != SMALLER);
      }
      else
      {
        // In this case the overlap is not valid.
        valid_intersection = false;
      }
    }
    is_first = false;

    if (valid_intersection)
      // Found an intersection to left_pt's right.
      break;

    // Discard the current intersection, which lies to left_pt's left.
    inter_list.pop_front();
  }

  // Insert the list of valid intersections into the map.
  inter_map[p_curve] = inter_list;

  // Return the first intersection object computed (may be empty).
  if (inter_list.empty())
    return CGAL::Object();
  else
    return (inter_list.front());
}

//-----------------------------------------------------------------------------
// Remove the next intersection of cv with the given halfedge from the map.
//
template<class Arrangement, class ZoneVisitor>
void Arrangement_zone_2<Arrangement,ZoneVisitor>::
_remove_next_intersection (Halfedge_handle he)
{
  // Get a pointer to the curve associated with the halfedge.
  const X_monotone_curve_2  *p_curve = &(he->curve());

  // Locate the intersections with this curve in the intersections map.
  Intersect_map_iterator     iter = inter_map.find (p_curve);

  CGAL_assertion (iter != inter_map.end());
  CGAL_assertion (! iter->second.empty());

  // Remove the first object in the list of intersections.
  iter->second.pop_front();
  return;
}

//-----------------------------------------------------------------------------
// Check if the given point lies completely to the left of the given egde.
//
template<class Arrangement, class ZoneVisitor>
bool Arrangement_zone_2<Arrangement,ZoneVisitor>::
_is_to_left_impl(const Point_2& p, Halfedge_handle he,
                 Arr_not_all_sides_oblivious_tag) const
{
  // Check the boundary conditions of the minimal end of the curve associated
  // with the given halfedge.
  const Arr_parameter_space   ps_x =
    m_geom_traits->parameter_space_in_x_2_object() (he->curve(), ARR_MIN_END);

  if (ps_x == ARR_LEFT_BOUNDARY)
    // The minimal end of the curve is to the left of any other point:
    return (false);

  const Arr_parameter_space   ps_y =
    m_geom_traits->parameter_space_in_y_2_object() (he->curve(), ARR_MIN_END);

  if (ps_y != ARR_INTERIOR) {
    // Check if p is to the left of the minimal curve-end:
    const Comparison_result   res =
      m_geom_traits->compare_x_point_curve_end_2_object() (p, he->curve(),
                                                           ARR_MIN_END);

    return ((res == SMALLER) || (res == EQUAL && ps_y == ARR_TOP_BOUNDARY));
  }

  // In case the minimal curve-end does not have boundary conditions, simply
  // compare p with the left endpoint of the curve.
  Vertex_const_handle   v_left = 
    (he->direction() == ARR_LEFT_TO_RIGHT) ? he->source() : he->target();

  return (m_geom_traits->compare_xy_2_object() (p, v_left->point()) == SMALLER);
}

//-----------------------------------------------------------------------------
// Check if the given point lies completely to the right of the given egde.
//
template<class Arrangement, class ZoneVisitor>
bool Arrangement_zone_2<Arrangement,ZoneVisitor>::
_is_to_right_impl(const Point_2& p, Halfedge_handle he,
                  Arr_not_all_sides_oblivious_tag) const
{
  // Check the boundary conditions of the maximal end of the curve associated
  // with the given halfedge.
  const Arr_parameter_space   ps_x =
    m_geom_traits->parameter_space_in_x_2_object() (he->curve(), ARR_MAX_END);

  if (ps_x == ARR_RIGHT_BOUNDARY)
    // The maximal end of the curve is to the right of any other point:
    return (false);

  const Arr_parameter_space   ps_y =
    m_geom_traits->parameter_space_in_y_2_object() (he->curve(), ARR_MAX_END);

  if (ps_y != ARR_INTERIOR) {
    // Check if p is to the right of the maximal curve-end:
    const Comparison_result   res =
      m_geom_traits->compare_x_point_curve_end_2_object() (p, he->curve(),
                                                           ARR_MAX_END);

    return ((res == LARGER) || (res == EQUAL && ps_y == ARR_BOTTOM_BOUNDARY));
  }

  // In case the maximal curve-end does not have boundary conditions, simply
  // compare p with the right endpoint of the curve.
  Vertex_const_handle   v_right = 
    (he->direction() == ARR_LEFT_TO_RIGHT) ? he->target() : he->source();

  return (m_geom_traits->compare_xy_2_object() (p, v_right->point()) == LARGER);
}

//-----------------------------------------------------------------------------
// Compute the (lexicographically) leftmost intersection of the query
// curve with the boundary of a given face in the arrangement.
//
template<class Arrangement, class ZoneVisitor>
void Arrangement_zone_2<Arrangement,ZoneVisitor>::
    _leftmost_intersection_with_face_boundary (Face_handle face,
                                               bool on_boundary)
{
  // Mark that we have not found any intersection (or overlap) yet.
  found_intersect = false;
  found_overlap = false;
  found_iso_vert = false;

  // Obtain some geometry-traits functors.
  typename Traits_adaptor_2::Compare_xy_2            compare_xy =
    m_geom_traits->compare_xy_2_object();
  typename Traits_adaptor_2::Is_in_x_range_2         is_in_x_range =
    m_geom_traits->is_in_x_range_2_object();
  typename Traits_adaptor_2::Construct_min_vertex_2  min_vertex =
    m_geom_traits->construct_min_vertex_2_object();
  typename Traits_adaptor_2::Compare_y_at_x_2        compare_y_at_x =
    m_geom_traits->compare_y_at_x_2_object();

  // Traverse the outer boundary of the face by going over all outer CCBs of
  // the face.
  typename Arrangement_2::Outer_ccb_iterator       occb_it;
  typename Arrangement_2::Ccb_halfedge_circulator  he_first;
  typename Arrangement_2::Ccb_halfedge_circulator  he_curr;

  CGAL::Object                 iobj;
  const Intersect_point_2     *int_p;
  const X_monotone_curve_2    *icv;
  Point_2                      ip;
  bool                         left_equals_curr_endpoint;
  bool                         intersection_on_right_boundary;
  bool                         leftmost_on_right_boundary = false;

  for (occb_it = face->outer_ccbs_begin();
       occb_it != face->outer_ccbs_end(); ++occb_it)
  {
    // Get circulators for the boundary of the current outer component.
    he_first = *occb_it;
    he_curr = he_first;

    do
    {
      // If this edge is fictitious, skip it.
      if (he_curr->is_fictitious())
      {
        ++he_curr;
        continue;
      }
      
      // If we have already found an intersection with the twin halfedge,
      // we do not have to compute intersections with the current halfedge.
      if (found_intersect && intersect_he == he_curr->twin())
      {
        ++he_curr;
        continue;
      }
      
      // If we already have an intersection point, compare it to the
      // endpoints of the curve associated with the current halfedge,
      // in order to filter unnecessary intersection computations.
      if (found_intersect && ! leftmost_on_right_boundary &&
          _is_to_left (intersect_p, he_curr))
      {
        // The current x-monotone curve lies entirely to the right of
        // ip_left, so its intersection with cv (if any) cannot lie to
        // the left of this point. We therefore do not need to compute
        // this intersection.
        ++he_curr;
        continue;
      }
      
      left_equals_curr_endpoint = false;
      if (on_boundary)
      {
        // Check if the left endpoint of the inserted curve (which is located
        // on the boundary of our face) equals one of the endpoints of the
        // current halfedge. If it equals the right endpoint of the current
        // halfedge, we can skip this edge, as there is no true overlap in
        // the x-range. Otherwise, we keep track of the fact that left_v is
        // the left end-vertex of the current halfedge.
        if (he_curr->target() == left_v)
        {
          left_equals_curr_endpoint = true;
          
          if (he_curr->direction() == ARR_LEFT_TO_RIGHT)
          {
            ++he_curr;
            continue;
          }
        }
        else if (he_curr->source() == left_v)
        {
          left_equals_curr_endpoint = true;
          
          if (he_curr->direction() == ARR_RIGHT_TO_LEFT)
          {
            ++he_curr;
            continue;
          }
        }
      }
    
      // Check whether the two curves overlap in their x-range (in order
      // to avoid unnecessary intersection computations).
      if (! left_equals_curr_endpoint &&
          ((! left_on_boundary && _is_to_right (left_pt, he_curr)) ||
           ! is_in_x_range (cv, he_curr->curve())))
      {
        // In case there is no overlap, the two x-monotone curves obviously
        // do not intersect.
        ++he_curr;
        continue;
      }

      // Compute the next intersection of cv and the current halfedge.
      iobj = _compute_next_intersection (he_curr,
                                         left_equals_curr_endpoint,
                                         intersection_on_right_boundary);

      if (! iobj.is_empty())
      {
        // We have found an intersection (either a simple point or an
        // overlapping x-monotone curve).
        int_p = object_cast<Intersect_point_2> (&iobj);
        if (int_p != NULL)
        {
          ip = int_p->first;

          // Found a simple intersection point. Check if it is the leftmost
          // intersection point so far.
          if (! found_intersect ||
              (! intersection_on_right_boundary &&
               (leftmost_on_right_boundary ||
                compare_xy (ip, intersect_p) == SMALLER)))
          {
            // Store the leftmost intersection point and the halfedge handle.
            intersect_p = ip;
            ip_mult = int_p->second;
            intersect_he = he_curr;
            found_overlap = false;
            leftmost_on_right_boundary = intersection_on_right_boundary;
          }
        }
        else
        {
          // We have located an overlapping curve. Assign ip as its left
          // endpoint.
          icv = object_cast<X_monotone_curve_2> (&iobj);
          CGAL_assertion (icv != NULL);

          ip = min_vertex (*icv);

          // Check if this endpoint it is the leftmost intersection point so
          // far.
          if (! found_intersect ||
              compare_xy (ip, intersect_p) == SMALLER)
          {
            // Store the leftmost intersection point and the halfedge handle.
            intersect_p = ip;
            ip_mult = 0;
            overlap_cv = *icv;
            intersect_he = he_curr;
            found_overlap = true;
          }
        }
        
        // Mark that we found an intersection.
        found_intersect = true;
      }
      
      // Move to the next edge along the outer boundary,
      ++he_curr;
    
    } while (he_curr != he_first);     // End loop on the current outer CCB.

  } // End: traversal of the outer CCBs of the face.

  // Traverse the inner boundary of the face by going over all inner CCBs
  // (the holes) of the face.
  typename Arrangement_2::Inner_ccb_iterator       iccb_it;

  for (iccb_it = face->inner_ccbs_begin();
       iccb_it != face->inner_ccbs_end(); ++iccb_it)
  {
    // Get circulators for the boundary of the current inner component.
    he_first = *iccb_it;
    he_curr = he_first;

    do
    {
      // If we have already found an intersection with the twin halfedge,
      // we do not have to compute intersections with the current halfedge.
      if (found_intersect && intersect_he == he_curr->twin())
      {
        ++he_curr;
        continue;
      }

      // If we already have an intersection point, compare it to the
      // endpoints of the curve associated with the current halfedge,
      //  in order to filter unnecessary intersection computations.
      if (found_intersect && ! leftmost_on_right_boundary &&
          _is_to_left (intersect_p, he_curr))
      {
        // The current x-monotone curve lies entirely to the right of
        // ip_left, so its intersection with cv (if any) cannot lie to
        // the left of this point. We therefore do not need to compute
        // this intersection.
        ++he_curr;
        continue;
      }

      left_equals_curr_endpoint = false;
      if (on_boundary)
      {
        // Check if the left endpoint of the inserted curve (which is located
        // on the boundary of our face) equals one of the endpoints of the
        // current halfedge. If it equals the right endpoint of the current
        // halfedge, we can skip this edge, as there is no true overlap in
        // the x-range. Otherwise, we keep track of the fact that left_v is
        // the left end-vertex of the current halfedge.
        if (he_curr->target() == left_v)
        {
          left_equals_curr_endpoint = true;

          if (he_curr->direction() == ARR_LEFT_TO_RIGHT)
          {
            ++he_curr;
            continue;
          }
        }
        else if (he_curr->source() == left_v)
        {
          left_equals_curr_endpoint = true;

          if (he_curr->direction() == ARR_RIGHT_TO_LEFT)
          {
            ++he_curr;
            continue;
          }
        }
      }

      // Check whether the two curves overlap in their x-range (in order
      // to avoid unnecessary intersection computations).
      if (! left_equals_curr_endpoint &&
          ((! left_on_boundary && _is_to_right (left_pt, he_curr)) ||
           ! is_in_x_range (cv, he_curr->curve())))
      {
        // In case there is no overlap, the two x-monotone curves obviously
        // do not intersect.
        ++he_curr;
        continue;
      }

      // Compute the next intersection of cv and the current halfedge.
      iobj = _compute_next_intersection (he_curr,
                                         left_equals_curr_endpoint,
                                         intersection_on_right_boundary);

      if (! iobj.is_empty())
      {
        // We have found an intersection (either a simple point or an
        // overlapping x-monotone curve).
        int_p = object_cast<Intersect_point_2> (&iobj);
        if (int_p != NULL)
        {
          ip = int_p->first;

          // Found a simple intersection point. Check if it is the leftmost
          // intersection point so far.
          if (! found_intersect ||
              (! intersection_on_right_boundary &&
               (leftmost_on_right_boundary ||
                compare_xy (ip, intersect_p) == SMALLER)))
          {
            // Store the leftmost intersection point and the halfedge
            // handle.
            intersect_p = ip;
            ip_mult = int_p->second;
            intersect_he = he_curr;
            found_overlap = false;
            leftmost_on_right_boundary = intersection_on_right_boundary;
          }
        }
        else
        {
          // We have located an overlapping curve. Assign ip as its left
          // endpoint.
          icv = object_cast<X_monotone_curve_2> (&iobj);
          CGAL_assertion (icv != NULL);

          ip = min_vertex (*icv);

          // Check if this endpoint it is the leftmost intersection point
          // so far.
          if (! found_intersect ||
              compare_xy (ip, intersect_p) == SMALLER)
          {
            // Store the leftmost intersection point and the halfedge
            // handle.
            intersect_p = ip;
            ip_mult = 0;
            overlap_cv = *icv;
            intersect_he = he_curr;
            found_overlap = true;
          }
        }

        // Mark that we found an intersection.
        found_intersect = true;
      }

      // Move to the next edge along the outer boundary,
      ++he_curr;

    } while (he_curr != he_first);      // End loop on the current inner CCB.

  } // End: traversal of the inner CCBs of the face.

  // Go over the boundary of the isolated vertices inside the face (if there
  // exist any), and check whether an isolated vertex lies on the curve.
  typename Arrangement_2::Isolated_vertex_iterator   iv_it;

  for (iv_it = face->isolated_vertices_begin();
       iv_it != face->isolated_vertices_end(); ++iv_it)
  {
    // If the isolated vertex is not in the x-range of our curve, disregard it.
    if (! is_in_x_range (cv, iv_it->point()))
      continue;

    // If we already have an intersection point, compare it to the current
    // isolated vertex, in order to filter unnecessary computations.
    if (found_intersect &&
        compare_xy (iv_it->point(), intersect_p) == LARGER)
    {
      continue;
    }

    // In case the isolated vertex lies on the curve, update the intersection
    // point accordingly.
    if (compare_y_at_x (iv_it->point(), cv) == EQUAL &&
        (! has_left_pt ||
         compare_xy (iv_it->point(), left_pt) == LARGER))
    {
      intersect_v = iv_it;
      intersect_p = intersect_v->point();
      ip_mult = 0;
      found_intersect = true;
      found_iso_vert = true;
    }

  } // End:: traversal of the isolated vertices inside the face.

  // Remove the next intersection associated with intersect_he, as we have
  // now reported it and do not want to encounter it again.
  if (found_intersect && !found_iso_vert)
    _remove_next_intersection (intersect_he);

  return;
}

//-----------------------------------------------------------------------------
// Compute the zone of an x-monotone curve in a given arrangement face.
// The left endpoint of the curve either lies in the face interior or on
// the boundary of the face.
//
template<class Arrangement, class ZoneVisitor>
bool Arrangement_zone_2<Arrangement,ZoneVisitor>::
_zone_in_face (Face_handle face, bool on_boundary)
{
  CGAL_precondition ((! on_boundary &&
                      ((left_v == invalid_v && left_he == invalid_he) ||
                      left_v->is_isolated())) ||
                     (on_boundary && left_he != invalid_he));

  // Find the first intersection of the curve with the face boundary.
  _leftmost_intersection_with_face_boundary (face, on_boundary);

  if (! found_intersect)
  {
    // Notify the visitor that the entire curve lies within the given face,
    // such that its right endpoint is not incident to any arrangement feature.
    visitor->found_subcurve (cv, face, left_v, left_he, invalid_v, invalid_he);

    // Inidicate that we are done with the zone-computation process.
    return (true);
  }

  // In this case found_intersect is true and intersect_he is the edge that
  // cv next intersects (or overlaps). If found_overlap is also true,
  // then overlap_cv is set and intersect_p is the left endpoint of the
  // overlapping subcurve. Otherwise, intersect_p is a simple intersection
  // point.
  // Alternatively, if found_iso_vert is true, then the next intersection point
  // intersect_p lies on the isolated vertex intersect_v.
  bool                  done = false;

  if (has_right_pt &&
      m_geom_traits->equal_2_object() (intersect_p, right_pt))
  {
    // If the intersection point is cv's right endpoint, the interior of cv
    // does not intersect any existing halfedge. In this case, we only have
    // to insert cv to the arrangement and we are done.
    sub_cv1 = cv;
    done = true;
  }
  else
  {
    // Split cv at the intersection point.
    m_geom_traits->split_2_object() (cv, intersect_p, sub_cv1, sub_cv2);

    // Set cv to be the remaining portion.
    has_left_pt = true;
    left_on_boundary = false;
    left_pt = intersect_p;
    cv = sub_cv2;
  }

  const X_monotone_curve_2  *p_orig_curve = NULL;

  if (! found_iso_vert)
  {
    // Check whether intersect_p coincides with one of the end-vertices of the
    // halfedge that cv intersects.
    if (! intersect_he->source()->is_at_open_boundary() &&
        m_geom_traits->equal_2_object() (intersect_p,
                                         intersect_he->source()->point()))
    {
      // We know that the right endpoint of sub_cv1 lies on the source vertex:
      right_v = intersect_he->source();
      right_he = invalid_he;
    }
    else if (! intersect_he->target()->is_at_open_boundary() &&
             m_geom_traits->equal_2_object() (intersect_p,
                                              intersect_he->target()->point()))
    {
      // We know that the right endpoint of sub_cv1 lies on the target vertex:
      right_v = intersect_he->target();
      right_he = invalid_he;
    }
    else
    {
      // The right endpoint of sub_cv1 lies on the interior of intersect_he:
      // Obtain the halfedge with the correct direction (which should be the
      // predecessor of sub_cv1 if we split the edge around this vertex).
      right_v = invalid_v;

      right_he = _direct_intersecting_edge_to_left (sub_cv1, intersect_he);
    }

    // Store the curve currently associated with the intersecting halfedge.
    p_orig_curve = &(intersect_he->curve());
  }
  else
  {
    // The right endpoint of the subcurve coincides with an isolated vertex:
    right_v = intersect_v;
    right_he = invalid_he;
  }

  // Notify the visitor that the left endpoint of the first subcurve is
  // located within the current face and both its endpoint are located
  // on its boundary.
  Visitor_result  visitor_res = visitor->found_subcurve (sub_cv1, face,
                                                         left_v, left_he,
                                                         right_v, right_he);

  // Check if we are done (either we have no remaining curve or if the
  // visitor has indicated we should end the process).
  if (done || visitor_res.second)
    return (true);

  // Move to the remaining portion of the curve, whose left endpoint is the
  // same as the right endpoint of sub_cv1. Note that we check if the visitor
  // has inserted the subcurve (in which case it should return a handle to
  // the resulting halfedge).
  Halfedge_handle  inserted_he = visitor_res.first;

  if (inserted_he != invalid_he)
  {
    if (right_v == invalid_v)
    {
      // If the right endpoint of the subcurve we have just detected was
      // not associated with an existing vertex, the inserted halfedge is
      // now targeted toward a newly created vertex that splits intersect_he
      // into two halfedges: (a) the next halfedge after inserted_he and (b)
      // the previous halfedge before inserted_he's twin.
      // The two halfedges (a) and (b) are now associated with the two
      // subcurves that result from splitting intersect_he->curve() at the
      // intersection point we have just detected, one extends to the left
      // and one to the right of this split point.
      const X_monotone_curve_2  *p_left_subcurve = NULL;
      const X_monotone_curve_2  *p_right_subcurve = NULL;

      if (inserted_he->next()->direction() == ARR_LEFT_TO_RIGHT)
      {
        // The next halfedge extends to the right of the split point:
        p_left_subcurve = &(inserted_he->twin()->prev()->curve());
        p_right_subcurve = &(inserted_he->next()->curve());
      }
      else
      {
        // The next halfedge extends to the left of the split point:
        p_right_subcurve = &(inserted_he->twin()->prev()->curve());
        p_left_subcurve = &(inserted_he->next()->curve());
      }

      // Associate the intersection list of the original curve with the
      // right subcurve, while we can associate an empty list with the
      // left subcurve, as we are now done with it.
      Intersect_map_iterator    iter = inter_map.find (p_orig_curve);
      Intersect_list            empty_inter_list;

      inter_map[p_right_subcurve] = iter->second;
      inter_map[p_left_subcurve] = empty_inter_list;

      // If necessary, erase the original curve from the intersection map.
      if (p_orig_curve != p_right_subcurve && p_orig_curve !=p_left_subcurve)
        inter_map.erase (p_orig_curve);
    }

    if (found_overlap && right_v == invalid_v)
    {
      // In case we have split the overlapping intersect_he, it now refers
      // to the wrong halfedge. the overlapping edge is either the successor
      // of the inserted halfedge or the predecessor of its twin, depending
      // on which one of these halfedges lies to the right of the split point.
      if (inserted_he->next()->direction() == ARR_LEFT_TO_RIGHT)
      {
        // The successor is directed to the right:
        intersect_he = inserted_he->next();
      }
      else
      {
        // The predecessor is directed to the left:
        CGAL_assertion (inserted_he->twin()->prev()->direction() ==
                        ARR_RIGHT_TO_LEFT);

        intersect_he = inserted_he->twin()->prev();
      }
    }

    // The visitor has created an edge that corresponds to sub_cv1 and instered
    // it into the arrangement. In this case, left_pt should be associated
    // with the target vertex of the new halfedge.
    CGAL_assertion
      (m_geom_traits->equal_2_object() (left_pt,
                                        inserted_he->target()->point()));

    left_v = inserted_he->target();

    // If right_he is known, it is possible to set left_he according to the
    // geometric information we have.
    if (right_he != invalid_he)
    {
      if ((ip_mult % 2) == 1)
      {
        // cv crosses right_he (which is now split into two), so the remaining
        // portion must be inserted after the next halfedge going clockwise
        // around left_v:
        //
        //              \   .                            .
        //               \ . remaining portion of cv     .
        //                x                              .
        //   inserted_he / \                             .
        //              /   \                            .
        left_he = inserted_he->next()->twin();
      }
      else if (ip_mult != 0)
      {
        // We have a tangency point. If right_he is directed from left to
        // right, we take the inserted halfedge to be left_he, otherwise
        // right_he itself becomes left_he:
        if (right_he->direction() == ARR_LEFT_TO_RIGHT)
          left_he = inserted_he;
        else
          left_he = right_he;
      }
      else
      {
        // Mutliplicity is unkown:
        left_he = invalid_he;
      }
    }
    else
    {
      // In case left_v used to be an isolated vertex, we know that the
      // inserted halfedge is its only incident halfedge and we can use it.
      // Otherwise, we do not know the identity of left_he.
      if (found_iso_vert)
        left_he = inserted_he;
      else
        left_he = invalid_he;
    }
  }
  else
  {
    // The visitor has not created a new edge. We proceed using the previously
    // computed arrangement features.
    left_v = right_v;

    if (right_he != invalid_he)
    {
      // In case cv crosses the interior of the right_he halfedge (the
      // multiplicity of the intersection is odd), we know that the ramaining
      // portion of the curve lies in the face incident to the twin halfedge.
      // If the multiplicity is known and is even, we stay with the same
      // halfedge.
      if ((ip_mult % 2) == 1)
        left_he = right_he->twin();
      else if (ip_mult != 0)
        left_he = right_he;
      else
        left_he = invalid_he;
    }
    else
    {
      left_he = invalid_he;
    }
  }

  // We are not done with the zone-computation process yet:
  return (false);
}

//-----------------------------------------------------------------------------
// Compute the zone of an overlapping subcurve overlap_cv of cv and the
// curve currently associated with intersect_he.
//
template<class Arrangement, class ZoneVisitor>
bool Arrangement_zone_2<Arrangement,ZoneVisitor>::_zone_in_overlap ()
{
  // Check if the right end of overlap_cv is bounded. If so, compute its
  // right endpoint.
  const bool           cv_has_right_pt =
    m_geom_traits->is_closed_2_object() (overlap_cv, ARR_MAX_END);

  Point_2              cv_right_pt;

  if (cv_has_right_pt)
    cv_right_pt = m_geom_traits->construct_max_vertex_2_object() (overlap_cv);

  // Get right end-vertex of the overlapping halfedge intersect_he. Also make
  // sure that the overlapping halfedge is always directed to the right.
  Vertex_handle   he_right_v;

  if (intersect_he->direction() == ARR_LEFT_TO_RIGHT)
  {
    he_right_v = intersect_he->target();
  }
  else
  {
    he_right_v = intersect_he->source();
    intersect_he = intersect_he->twin();
  }

  // Compare the two right endpoints. Note that overlap_cv cannot extend to
  // the right longer than the halfedge it overlaps. Thus, if the curve is
  // not bounded, the right vertex of intersect_he must lie on open boundary as
  // well.
  if (! cv_has_right_pt)
  {
    CGAL_assertion_code
      (const Arr_parameter_space  cv_ps_x = 
       m_geom_traits->parameter_space_in_x_2_object() (overlap_cv, ARR_MAX_END);
       const Arr_parameter_space  cv_ps_y =
       m_geom_traits->parameter_space_in_y_2_object() (overlap_cv, ARR_MAX_END);
       );
    CGAL_assertion (he_right_v->parameter_space_in_x() == cv_ps_x &&
                    he_right_v->parameter_space_in_y() == cv_ps_y);

    right_v = he_right_v;
  }
  else
  {
    // In this case overlap_cv has a finite right endpoint. In this case,
    // if the right vertex of intersect_he is associated with a finite point,
    // we check whether it is equal to cv_right_pt. Otherwise, we know that
    // intersect_he extends to the the right of overlap_cv, and there is no
    // vertex currently associated with overlap_cv's right endpoint.
    if (! he_right_v->is_at_open_boundary() &&
        m_geom_traits->equal_2_object() (cv_right_pt, he_right_v->point()))
    {
      // The overlap is with the entire halfedge. In this case we set the
      // right end-vertex of the overlapping zone.
      right_v = he_right_v;
    }
    else
    {
      // In this case intersect_he overlaps just a portion of prev_he.
      // The right end-vertex of the overlapping zone is not known.
      right_v = invalid_v;
    }
  }

  // Store the curve currently associated with the overlapping halfedge.
  const X_monotone_curve_2  *p_orig_curve = &(intersect_he->curve());

  // Notify the visitor on the overlapping zone.
  Visitor_result  visitor_res = visitor->found_overlap (overlap_cv,
                                                        intersect_he,
                                                        left_v, right_v);

  // If the visitor has indicated we should halt the process, or it the right
  // endpoint of the overlapping curve is the right endpoint of cv then we are
  // done (or both extend to an open boundary).
  if (visitor_res.second ||
      (cv_has_right_pt && has_right_pt &&
       m_geom_traits->equal_2_object() (cv_right_pt, right_pt)) ||
      (! cv_has_right_pt && ! has_right_pt))
  {
    return (true);
  }

  // Erase the original curve from the intersection map, so we will have to
  // recompute intersections with it in the future.
  inter_map.erase (p_orig_curve);

  // Mark that we have dealt with the overlap.
  found_overlap = false;

  // Split cv at right endpoint of the overlapping curve.
  m_geom_traits->split_2_object() (cv, cv_right_pt, sub_cv1, sub_cv2);

  // Set cv to be the remaining portion.
  has_left_pt = true;
  left_on_boundary = false;
  left_pt = cv_right_pt;
  cv = sub_cv2;

  // Move to the remaining portion of the curve, whose left endpoint is the
  // same as the right endpoint of the overlapping curve. Note that we check
  // if the visitor has inserted the subcurve (in which case it should return
  // a handle to the resulting halfedge).
  Halfedge_handle  updated_he = visitor_res.first;

  if (updated_he != invalid_he)
  {
    // In this case, left_pt should be associated with the target vertex of
    // the updated halfedge.
    CGAL_assertion
      (m_geom_traits->equal_2_object() (left_pt, updated_he->target()->point()));
 
    left_v = updated_he->target();
  }
  else
  {
    left_v = right_v;
  }

  left_he = invalid_he;

  // We are not done with the zone-computation process yet:
  return (false);
}

} //namespace CGAL

#endif
