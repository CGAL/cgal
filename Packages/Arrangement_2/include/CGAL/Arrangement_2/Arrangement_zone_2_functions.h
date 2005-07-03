// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Ron Wein          <wein@post.tau.ac.il>
//                 (based on old version by Eyal Flato)

#ifndef CGAL_ARRANGEMENT_ZONE_2_FUNCTIONS_H
#define CGAL_ARRANGEMENT_ZONE_2_FUNCTIONS_H

/*! \file
 * Member-function definitions for the Arrangement_zone_2 class.
 */

CGAL_BEGIN_NAMESPACE

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
  typename Arrangement_2::Vertex_const_handle    vh;
  typename Arrangement_2::Halfedge_const_handle  hh;
  typename Arrangement_2::Face_const_handle      fh;

  if (assign (vh, obj))
  {
    // The left endpoint coincides with an existing vertex:
    left_v = arr.non_const_handle (vh);
  }
  else if (assign (hh, obj))
  {
    // Obtain the right halfedge from the halfedge-pair containing left_pt
    // in their interior.
    left_he = _direct_intersecting_edge_to_right (cv, left_pt,
                                                  arr.non_const_handle (hh));

    // Handle overlaps.
    if (found_overlap)
    {
      // In this case cv overlaps the curve associated with intersect_he.
      // Compute the overlapping subcurve to the right of curr_v.
      obj = _compute_next_intersection (intersect_he);

      bool  assign_success = assign (overlap_cv, obj);

      CGAL_assertion (assign_success);
      if (assign_success)
      {
        // Remove the overlap from the map.
        _remove_next_intersection (intersect_he);

        // Compute the overlap zone and continue to the end of the loop.
        done = _zone_in_overlap ();
      }
    }
  }
  else
  {
    // The left endpoint lies inside a face.
    bool     assign_success = assign (fh, obj);

    CGAL_assertion_msg(assign_success,
		       "Invalid object returned by the point-location query.");

    if (! assign_success)
      return;

    // Compute the zone of the curve at the interior of the face.
    done = _zone_in_face (arr.non_const_handle(fh),
                          false);      // left_pt is not on the face boundary.

    // In case we have just discovered an overlap, compute the overlapping
    // xone as well.
    if (! done && found_overlap)
      done = _zone_in_overlap ();
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
          found_overlap = _find_prev_around_vertex (left_v,
                                                    left_he);
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
        left_he = _direct_intersecting_edge_to_right (cv, left_pt,
                                                      right_he);
      }

      if (found_overlap)
      {
        // In this case cv overlaps the curve associated with intersect_he.
        // Compute the overlapping subcurve to the right of curr_v.
        obj = _compute_next_intersection (intersect_he);

        bool  assign_success = assign (overlap_cv, obj);

        CGAL_assertion (assign_success);
        if (assign_success)
        {
          // Remove the overlap from the map.
          _remove_next_intersection (intersect_he);

          // Compute the overlap zone and continue to the end of the loop.
          done = _zone_in_overlap ();
          continue;
        }
      }
    }

    if (left_v == invalid_v || ! left_v->is_isolated())
    {
      // At this point we can compute the zone of cv starting from the left_he
      // inside its incident face.
      done = _zone_in_face (left_he->face(),
                            true);      // left_pt is on the face boundary.
    }
    else
    {
      // Compute the zone of cv starting from the face that contains the
      // isolated vertex left_v.
      done = _zone_in_face (arr.incident_face(left_v),
                            false);     // left_pt is not on the face boundary.
    }

    // In case we have just discovered an overlap, compute the overlapping
    // zone as well.
    if (! done && found_overlap)
    {
      // Remove the overlap from the map.
      _remove_next_intersection (intersect_he);

      done = _zone_in_overlap();
    }
  }

  // Clear the intersections map and the list of halfedges that originated
  // from cv.
  inter_map.clear();
  invalid_cvs.clear();

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
    he = he_curr->handle();

    // Note that cv extends to the right of v. In case the single
    // halfedge also extends to the right of v (its source is to
    // the right), check if an overlap occurs.
    if ((traits->compare_xy_2_object()(v->point(),
                                       he->source()->point()) == SMALLER) &&
        (traits->compare_y_at_x_right_2_object() (he->curve(), cv,
                                                  v->point()) == EQUAL))
    {
      // Mark that an overlap occurs:
      intersect_he = he_curr->handle();
      return (true);
    }

    // We have no overlap:
    return (false);
  }

  // Find the face containing cv around the vertex.
  typename Traits_wrapper_2::Is_between_cw_2  is_between_cw =
                                             traits->is_between_cw_2_object();

  do
  {
    // Check if it is possible to insert cv in between the current curve
    // and the next curve, going in a clockwise direction around v.
    is_between = is_between_cw (cv,
                                he_curr->curve(), he_next->curve(),
                                v->point(),
                                cv_equals_curr, cv_equals_next);

    // Check the case of overlaps:
    if (cv_equals_curr)
    {
      // cv overlaps with the curve of he_curr:
      intersect_he = he_curr->handle();
      return (true);
    }
    else if (cv_equals_next)
    {
      // cv overlaps with the curve of he_next:
      intersect_he = he_next->handle();
      return (true);
    }

    if (is_between)
    {
      // We can conclude that cv should be placed between he_curr and
      // he_next (in a clockwise order), and no overlap occurs.
      he = he_curr->handle();
      return (false);
    }

    // Proceed to the next halfedges around the vertex.
    ++he_curr;
    ++he_next;

  } while (he_curr != he_first);

  // We should never reach here:
  CGAL_assertion (false);
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
  CGAL_assertion (traits->compare_y_at_x_2_object()
                    (cv_left_pt, query_he->curve()) == EQUAL);

  // Check whether the given halfedge is directed to the right.
  const bool               query_he_directed_right =
    (traits->compare_xy_2_object() (query_he->source()->point(),
                                    query_he->target()->point()) == SMALLER);

  // Check whether the curve lies above of below the edge immediately to
  // the right of its left endpoint.
  const Comparison_result  pos_res =
    traits->compare_y_at_x_right_2_object() (cv_ins, query_he->curve(),
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
Arrangement_zone_2<Arrangement,ZoneVisitor>::_direct_intersecting_edge_to_left
    (const X_monotone_curve_2& cv_ins,
     Halfedge_handle query_he)
{
  // Make sure that the right endpoint of cv_ins lies on query_he.
  CGAL_assertion (traits->compare_y_at_x_2_object()
                  (traits->construct_max_vertex_2_object() (cv_ins),
                   query_he->curve()) == EQUAL);

  // Check whether the given halfedge is directed to the right.
  const bool               query_he_directed_right =
    (traits->compare_xy_2_object() (query_he->source()->point(),
                                    query_he->target()->point()) == SMALLER);

  // Check whether the curve lies above of below the edge (we use the curve
  // position predicate, as we know they cruves do not overlap and intersect
  // only at the split point).
  const Comparison_result  pos_res =
      traits->compare_y_position_2_object() (cv_ins, query_he->curve());

  if (pos_res == SMALLER)
  {
    // If cv_ins lies below the curve associated with query_he, we should
    // take the haldedge directed from right to left, so if query_he is
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
    // take the haldedge directed from left to right, so if query_he is
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
Object Arrangement_zone_2<Arrangement,ZoneVisitor>::_compute_next_intersection
    (Halfedge_handle he)
{
  // Get a pointer to the curve associated with the halfedge.
  const X_monotone_curve_2  *p_curve = &(he->curve());

  // Check if the curve belongs to the set of invalid curves, with whom the
  // previously computed intersections are no longer valid.
  Curves_set_iterator        cv_iter = invalid_cvs.find (p_curve);
  bool                       invalid_intersections_with_curve = false;

  if (cv_iter != invalid_cvs.end())
  {
    invalid_intersections_with_curve = true;
  }

  // Try to locate the intersections with this curve in the intersections map.
  Intersect_map_iterator     iter;

  if (! invalid_intersections_with_curve)
    iter = inter_map.find (p_curve);
  else
    iter = inter_map.end();

  if (iter != inter_map.end())
  {
    // The intersections with the curve have already been computed.
    // Retrieve the intersections list from the map.
    Intersect_list&          inter_list = iter->second;

    if (inter_list.empty())
      return Object();

    // Locate the first intersection that lies to the right of left_pt.
    Intersect_point_2        ip;
    bool                     valid_intersection;

    while (! inter_list.empty())
    {
      // Compare that current object with left_pt.
      if (assign (ip, inter_list.front()))
      {
        valid_intersection =
          (traits->compare_xy_2_object() (ip.first, left_pt) == LARGER);
      }
      else
      {
        X_monotone_curve_2   icv;

        assign (icv, inter_list.front());
        valid_intersection = (traits->compare_xy_2_object()
                              (traits->construct_min_vertex_2_object()(icv),
                               left_pt) != SMALLER);
      }

      if (valid_intersection)
        // Found an intersection to left_pt's right.
        return (inter_list.front());

      // Discard the current intersection, which lies to left_pt's left.
      inter_list.pop_front();
    }

    // If we reached here, the list of intersections is empty:
    return Object();
  }

  // The intersections with the curve have not been computed yet, so we
  // have to compute them now.
  Intersect_list           inter_list;

  traits->intersect_2_object() (cv, he->curve(),
                                std::back_inserter(inter_list));

  // Discard all intersection lying to the left of left_pt.
  Intersect_point_2        ip;
  bool                     valid_intersection;

  while (! inter_list.empty())
  {
    // Compare that current object with left_pt.
    if (assign (ip, inter_list.front()))
    {
      valid_intersection =
        (traits->compare_xy_2_object() (ip.first, left_pt) == LARGER);
    }
    else
    {
      X_monotone_curve_2   icv;

      assign (icv, inter_list.front());
      valid_intersection = (traits->compare_xy_2_object()
                            (traits->construct_min_vertex_2_object() (icv),
                             left_pt) != SMALLER);
    }

    if (valid_intersection)
      // Found an intersection to left_pt's right.
      break;

    // Discard the current intersection, which lies to left_pt's left.
    inter_list.pop_front();
  }

  // Insert the list of valid intersections into the map.
  inter_map[p_curve] = inter_list;

  // In case the curve used to belong to the set of invalid curves, remove it
  // from there, as we have now updated the intersections list.
  if (cv_iter != invalid_cvs.end())
    invalid_cvs.erase (cv_iter);

  // Return the first intersection object computed (may be empty).
  if (inter_list.empty())
    return Object();
  else
    return (inter_list.front());
}

//-----------------------------------------------------------------------------
// Remove the next intersection of cv with the given halfedge from the map.
//
template<class Arrangement, class ZoneVisitor>
void Arrangement_zone_2<Arrangement,ZoneVisitor>::_remove_next_intersection
    (Halfedge_handle he)
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
// Compute the (lexicographically) leftmost intersection of the query
// curve with the boundary of a given face in the arrangement.
//
template<class Arrangement, class ZoneVisitor>
void Arrangement_zone_2<Arrangement,ZoneVisitor>::
    _leftmost_intersection_with_face_boundary (Face_handle face)
{
  // Mark that we have not found any intersection (or overlap) yet.
  found_intersect = false;
  found_overlap = false;
  found_iso_vert = false;

  // Go over the outer boundary of the face (if one exists), and try to
  // locate intersections of cv with the edges along the boundary.
  typename Traits_wrapper_2::Compare_xy_2            compare_xy =
                                      traits->compare_xy_2_object();
  typename Traits_wrapper_2::Is_in_x_range_2         is_in_x_range =
                                      traits->is_in_x_range_2_object();
  typename Traits_wrapper_2::Construct_min_vertex_2  min_vertex =
                                      traits->construct_min_vertex_2_object();
  typename Traits_wrapper_2::Compare_y_at_x_2        compare_y_at_x =
                                      traits->compare_y_at_x_2_object();

  typename Arrangement_2::Ccb_halfedge_circulator  he_first;
  typename Arrangement_2::Ccb_halfedge_circulator  he_curr;

  Object                  obj;
  Intersect_point_2       int_p;
  Point_2                 ip;

  if (! face->is_unbounded())
  {
    // Get circulators for the outer boundary of the face.
    he_first = face->outer_ccb();
    he_curr = he_first;

    do
    {
      // If we already have an intersection point, compare it to the
      // endpoints of the curve associated with the current halfedge,
      //  in order to filter unnecessary intersection computations.
      if (found_intersect &&
          compare_xy (he_curr->source()->point(), intersect_p) == LARGER &&
          compare_xy (he_curr->target()->point(), intersect_p) == LARGER)
      {
        // The current x-monotone curve lies entirely to the right of
        // ip_left, so its intersection with cv (if any) cannot lie to
        // the left of this point. We therefore do not need to compute
        // this intersection.
        ++he_curr;
        continue;
      }

      // Check whether the two curves overlap in their x-range (in order
      // to avoid unnecessary intersection computations).
      if (! is_in_x_range (cv, he_curr->curve()))
      {
        // In case there is no overlap, the two x-monotone curve obviously
        // do not intersect.
        ++he_curr;
        continue;
      }

      // Compute the next intersection of cv and the current halfedge.
      obj = _compute_next_intersection (he_curr);

      if (! obj.is_empty())
      {
        // We have found an intersection (either a simple point or an
        // overlapping x-monotone curve).
        if (assign (int_p, obj))
        {
          ip = int_p.first;

          // Found a simple intersection point. Check if it is the leftmost
          // intersection point so far.
          if (! found_intersect ||
              compare_xy (ip, intersect_p) == SMALLER)
          {
            // Store the leftmost intersection point and the halfedge handle.
            intersect_p = ip;
            ip_mult = int_p.second;
            intersect_he = he_curr->handle();
            found_overlap = false;
          }
        }
        else
        {
          // We have located an overlapping curve. Assign ip as its left
          // endpoint.
          X_monotone_curve_2   icv;

          assign (icv, obj);
          ip = min_vertex (icv);

          // Check if this endpoint it is the leftmost intersection point so
          // far.
          if (! found_intersect ||
              compare_xy (ip, intersect_p) == SMALLER)
          {
            // Store the leftmost intersection point and the halfedge handle.
            intersect_p = ip;
            ip_mult = 0;
            overlap_cv = icv;
            intersect_he = he_curr->handle();
            found_overlap = true;
          }
        }

        // Mark that we found an intersection.
        found_intersect = true;
      }

      // Move to the next edge along the outer boundary,
      ++he_curr;

    } while (he_curr != he_first);

  } // End: if (! face->is_unbounded())

  // Go over the boundary of the holes inside the face (if there exist any),
  // and try to locate intersections of cv with the edges along the boundary
  // of each hole.
  typename Arrangement_2::Holes_iterator   holes_it;

  for (holes_it = face->holes_begin();
       holes_it != face->holes_end(); ++holes_it)
  {
    // Get circulators for the boundary of the current hole.
    he_first = *holes_it;
    he_curr = he_first;

    do
    {
      // If we already have an intersection point, compare it to the
      // endpoints of the curve associated with the current halfedge,
      //  in order to filter unnecessary intersection computations.
      if (found_intersect &&
          compare_xy (he_curr->source()->point(), intersect_p) == LARGER &&
          compare_xy (he_curr->target()->point(), intersect_p) == LARGER)
      {
        // The current x-monotone curve lies entirely to the right of
        // ip_left, so its intersection with cv (if any) cannot lie to
        // the left of this point. We therefore do not need to compute
        // this intersection.
        ++he_curr;
        continue;
      }

      // Check whether the two curves overlap in their x-range (in order
      // to avoid unnecessary intersection computations).
      if (! is_in_x_range (cv, he_curr->curve()))
      {
        // In case there is no overlap, the two x-monotone curve obviously
        // do not intersect.
        ++he_curr;
        continue;
      }

      // Compute the next intersection of cv and the current halfedge.
      obj = _compute_next_intersection (he_curr);

      if (! obj.is_empty())
      {
        // We have found an intersection (either a simple point or an
        // overlapping x-monotone curve).
        if (assign (int_p, obj))
        {
          ip = int_p.first;

          // Found a simple intersection point. Check if it is the leftmost
          // intersection point so far.
          if (! found_intersect ||
              compare_xy (ip, intersect_p) == SMALLER)
          {
            // Store the leftmost intersection point and the halfedge
            // handle.
            intersect_p = ip;
            ip_mult = int_p.second;
            intersect_he = he_curr->handle();
            found_overlap = false;
          }
        }
        else
        {
          // We have located an overlapping curve. Assign ip as its left
          // endpoint.
          X_monotone_curve_2   icv;

          assign (icv, obj);
          ip = min_vertex (icv);

          // Check if this endpoint it is the leftmost intersection point
          // so far.
          if (! found_intersect ||
              compare_xy (ip, intersect_p) == SMALLER)
          {
            // Store the leftmost intersection point and the halfedge
            // handle.
            intersect_p = ip;
            ip_mult = 0;
            overlap_cv = icv;
            intersect_he = he_curr->handle();
            found_overlap = true;
          }
        }

        // Mark that we found an intersection.
        found_intersect = true;
      }

      // Move to the next edge along the outer boundary,
      ++he_curr;

    } while (he_curr != he_first);

  } // End: traversal of the holes inside the face.

  // Go over the boundary of the isolated vertices inside the face (if there
  // exist any), and check whether an isolated vertex lies on the curve.
  typename Arrangement_2::Isolated_vertices_iterator   iso_verts_it;

  for (iso_verts_it = face->isolated_vertices_begin();
       iso_verts_it != face->isolated_vertices_end(); ++iso_verts_it)
  {
    // If the isolated vertex is not in the x-range of our curve, disregard it.
    if (! is_in_x_range (cv, iso_verts_it->point()))
      continue;

    // If we already have an intersection point, compare it to the current
    // isolated vertex, in order to filter unnecessary computations.
    if (found_intersect &&
        compare_xy (iso_verts_it->point(), intersect_p) == LARGER)
    {
      continue;
    }

    // In case the isolated vertex lies on the curve, update the intersection
    // point accordingly.
    if (compare_y_at_x (iso_verts_it->point(), cv) == EQUAL &&
        compare_xy (iso_verts_it->point(), left_pt) == LARGER)
    {
      intersect_v = iso_verts_it->handle();
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
bool Arrangement_zone_2<Arrangement,ZoneVisitor>::_zone_in_face
    (Face_handle face,
     bool on_boundary)
{
  CGAL_precondition ((! on_boundary &&
                      ((left_v == invalid_v && left_he == invalid_he) ||
                      left_v->is_isolated())) ||
                     (on_boundary && left_he != invalid_he));

  // Find the first intersection of the curve with the face boundary.
  _leftmost_intersection_with_face_boundary (face);

  if (! found_intersect)
  {
    // Notify the visitor that the entire curve lies within the given face,
    // such that its right endpoint is not incident to any arrangement feature.
    visitor->found_subcurve (cv,
                             face,
                             left_v, left_he,
                             invalid_v, invalid_he);

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

  if (traits->equal_2_object() (intersect_p, right_pt))
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
    traits->split_2_object() (cv,
                              intersect_p,
                              sub_cv1, sub_cv2);

    // Set cv to be the remaining portion.
    left_pt = intersect_p;
    cv = sub_cv2;
  }

  if (! found_iso_vert)
  {
    // Check whether intersect_p coincides with one of the end-vertices of the
    // halfedge that cv intersects.
    if (traits->equal_2_object() (intersect_p,
                                  intersect_he->source()->point()))
    {
      // We know that the right endpoint of sub_cv1 lies on the source vertex:
      right_v = intersect_he->source();
      right_he = invalid_he;
    }
    else if (traits->equal_2_object() (intersect_p,
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
  }
  else
  {
    // The right endpoint of the subcurve coincides with an isolated vertex:
    right_v = intersect_v;
    right_he = invalid_he;
  }

  // Store the curve currently associated with the intersecting halfedge.
  const X_monotone_curve_2  *p_old_curve = &(intersect_he->curve());

  // Notify the visitor that the left endpoint of the first subcurve is
  // located within the current face and both its endpoint are located
  // on its boundary.
  Visitor_result  visitor_res = visitor->found_subcurve (sub_cv1,
                                                         face,
                                                         left_v, left_he,
                                                         right_v, right_he);

  // Check if we are done (either we have no remaining curve or if the
  // visitor has indicated we should end the process).
  if (done || visitor_res.second)
    return (true);

  // Mark that the curve that used to be associated with the intersecting
  // halfedge becomes invalid now.
  invalid_cvs.insert (p_old_curve);

  // Move to the remaining portion of the curve, whose left endpoint is the
  // same as the right endpoint of sub_cv1. Note that we check if the visitor
  // has inserted the subcurve (in which case it should return a handle to
  // the resulting halfedge).
  Halfedge_handle  inserted_he = visitor_res.first;

  if (inserted_he != invalid_he)
  {
    // The visitor has created an edge that corresponds to sub_cv1 and instered
    // it into the arrangement. In this case, left_pt should be associated
    // with the target vertex of the new halfedge.
    CGAL_assertion (traits->equal_2_object() (left_pt,
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
        if (traits->compare_xy_2_object()
            (right_he->source()->point(),
             right_he->target()->point()) == SMALLER)
        {
          left_he = inserted_he;
        }
        else
        {
          left_he = right_he;
        }
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
  // Get the right endpoints of overlap_cv and the right end-vertex of
  // the overlapping halfedge intersect_he. Also make sure that the overlapping
  // halfedge is always directed to the right.
  Point_2         cv_right_pt =
    traits->construct_max_vertex_2_object() (overlap_cv);
  Vertex_handle   he_right_v;

  if (traits->compare_xy_2_object() (intersect_he->source()->point(),
                                     intersect_he->target()->point()) == SMALLER)
  {
    he_right_v = intersect_he->target();
  }
  else
  {
    he_right_v = intersect_he->source();
    intersect_he = intersect_he->twin();
  }

  // Compare the two right endpoints. Note that overlap_cv cannot extend to
  // the right longer than the halfedge it overlaps.
  const Point_2&         he_right_pt = he_right_v->point();

  Comparison_result      res = traits->compare_xy_2_object() (cv_right_pt,
                                                              he_right_pt);

  CGAL_assertion (res != LARGER);

  if (res == EQUAL)
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

  // Store the curve currently associated with the overlapping halfedge.
  const X_monotone_curve_2  *p_old_curve = &(intersect_he->curve());

  // Notify the visitor on the overlapping zone.
  Visitor_result  visitor_res = visitor->found_overlap (overlap_cv,
                                                        intersect_he,
                                                        left_v, right_v);

  // If the visitor has indicated we should halt the process, or it the right
  // endpoint of the overlapping curve is the right endpoint of cv then we are
  // done.
  if (visitor_res.second ||
      traits->equal_2_object() (cv_right_pt, right_pt))
  {
    return (true);
  }

  // Mark that the curve that used to be associated with the overlapping
  // halfedge becomes invalid now.
  invalid_cvs.insert (p_old_curve);

  // Mark that we have dealt with the overlap.
  found_overlap = false;

  // Split cv at right endpoint of the overlapping curve.
  traits->split_2_object() (cv,
                            cv_right_pt,
                            sub_cv1, sub_cv2);

  // Set cv to be the remaining portion.
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
    CGAL_assertion (traits->equal_2_object() (left_pt,
                                              updated_he->target()->point()));
 
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

CGAL_END_NAMESPACE

#endif
