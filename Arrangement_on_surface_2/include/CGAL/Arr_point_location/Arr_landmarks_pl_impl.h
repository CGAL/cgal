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
// $URL$
// $Id$
//
// Author(s)     : Idit Haran   <haranidi@post.tau.ac.il>
//                 Ron Wein     <wein@post.tau.ac.il>
#ifndef CGAL_ARR_LANDMARKS_PL_IMPL_H
#define CGAL_ARR_LANDMARKS_PL_IMPL_H

/*! \file
 * Member-function definitions for the 
 * Arr_landmarks_point_location<Arrangement, Generator> class.
 */

CGAL_BEGIN_NAMESPACE

//-----------------------------------------------------------------------------
// Locate the arrangement feature containing the given point.
//
template <class Arr, class Gen>
Object Arr_landmarks_point_location<Arr, Gen>::locate
    (const Point_2& p) const
{
  // If the arrangement is empty, return its initial (empty and 
  // non-fictitious) face.
  if (p_arr->number_of_vertices() == 0)
  {
    CGAL_assertion (p_arr->number_of_faces() == 1);

    Face_const_handle      fh = p_arr->faces_begin();

    return (CGAL::make_object (fh));
  }

  // Use the generator and to find the closest landmark to the query point.
  Object         lm_location_obj; 
  const Point_2& landmark_point = lm_gen->closest_landmark (p, lm_location_obj);

  if (landmark_point==p)
  {
    //the easy case, the query point and the landmark point
    //are the same so return the landmark location
    return lm_location_obj;
  }

  // Walk from the nearest_vertex to the point p, using walk algorithm, 
  // and find the location of the query point p. Note that the set fo edges
  // we have crossed so far is initially empty.
  Halfedge_set                   crossed_edges;
  Object                         out_obj;
  
  // Locate the arrangement feature that contains the landmark.
  const Vertex_const_handle     *vh;
  const Halfedge_const_handle   *hh;
  const Face_const_handle       *fh;

  CGAL_assertion_msg (! lm_location_obj.is_empty(),
                      "lm_location_obj is empty.");

  if ((vh = object_cast<Vertex_const_handle>(&lm_location_obj)) != NULL)
  {
    out_obj = _walk_from_vertex (*vh, p, crossed_edges);
  }
  else if ((fh = object_cast<Face_const_handle>(&lm_location_obj)) != NULL)
  {
    out_obj = _walk_from_face (*fh, landmark_point, p, crossed_edges);
  }
  else if ((hh = object_cast<Halfedge_const_handle>(&lm_location_obj)) != NULL)
  {
    out_obj = _walk_from_edge (*hh, landmark_point, p, crossed_edges);
  }
  else 
  {
    CGAL_assertion_msg (! lm_location_obj.is_empty(),
                        "lm_location_obj of an unknown type.");
    return (out_obj);
  }

  if ((fh = object_cast<Face_const_handle>(&out_obj)) != NULL)
  {
    // If we reached here, we did not locate the query point in any of the
    // holes inside the current face, so we conclude it is contained in this
    // face.
    // However, we first have to check whether the query point coincides with
    // any of the isolated vertices contained inside this face.
    Isolated_vertex_const_iterator      iso_verts_it;
    typename Traits_adaptor_2::Equal_2  equal = m_traits->equal_2_object();

    for (iso_verts_it = (*fh)->isolated_vertices_begin();
         iso_verts_it != (*fh)->isolated_vertices_end(); ++iso_verts_it)
    {
      if (equal (p, iso_verts_it->point()))
      {
        Vertex_const_handle  ivh = iso_verts_it;
        return (CGAL::make_object (ivh));
      }
    }
  }

  return (out_obj);
}

//-----------------------------------------------------------------------------
// Walk from a given vertex to the query point.
//
template <class Arr, class Gen>
Object Arr_landmarks_point_location<Arr, Gen>::_walk_from_vertex
    (Vertex_const_handle nearest_vertex,
     const Point_2& p,
     Halfedge_set& crossed_edges) const
{ 
  Vertex_const_handle       vh = nearest_vertex;

  CGAL_assertion_msg (! vh->is_at_infinity(),
                      "_walk_from_vertex() from a vertex at infinity.");

  // Check if the qurey point p conincides with the vertex.
  if (m_traits->equal_2_object() (vh->point(), p))
    return (CGAL::make_object (vh));

  // In case of an isolated vertex, walk to from the face that contains
  // it toward the query point.
  if (vh->is_isolated())
  {
    Face_const_handle  fh = vh->face();
    return (_walk_from_face (fh, vh->point(), p, crossed_edges));
  }

  // Locate the face around the vertex that contains the curve connecting
  // the vertex and the query point.
  bool                      new_vertex = false;
  Object                    obj;
  const Face_const_handle  *p_fh;

  do
  {
    new_vertex = false;
    obj = _find_face_around_vertex (vh, p, new_vertex);

    if (new_vertex)
    {
      // We found a vertex closer to p: Continue using this vertex.
      vh = object_cast<Vertex_const_handle> (obj);
    }
    else if (object_cast<Halfedge_const_handle>(&obj) != NULL ||
             object_cast<Vertex_const_handle>(&obj) != NULL)
    {
      // If p is located on an edge or on a vertex, return the object
      // that wraps this arrangement feature.
      return (obj);
    }
    else if ((p_fh = object_cast<Face_const_handle>(&obj)) != NULL)
    {
      // Walk to p from the face we have located:
      return (_walk_from_face (*p_fh, vh->point(), p, crossed_edges));
    }
    else
    {
      CGAL_assertion_msg
          (false, "_find_face_around_vertex() returned an unknown object.");
    }

  } while (new_vertex);

  // We should never reach here:
  CGAL_error();
  return Object();
}

//-----------------------------------------------------------------------------
// Locate an edge around a given vertex that is the predecessor of the curve
// connecting the vertex to the query point in a clockwise order.
// 
template <class Arr, class Gen>
Object Arr_landmarks_point_location<Arr, Gen>::_find_face_around_vertex
    (Vertex_const_handle vh,
     const Point_2& p, 
     bool& new_vertex) const
{
  new_vertex = false;

  // Create an x-monotone curve connecting the point associated with the
  // vertex vp and the query point p.
  const Point_2&      vp = vh->point();
  X_monotone_curve_2  seg =
    m_traits->construct_x_monotone_curve_2_object()(vp, p);
  const bool          seg_dir_right =
    (m_traits->compare_xy_2_object()(vp, p) == SMALLER);

  // Get the first incident halfedge around v and the next halfedge.
  Halfedge_around_vertex_const_circulator  first = vh->incident_halfedges();
  Halfedge_around_vertex_const_circulator  curr, next;
  bool                                     equal_curr = false;

  next = curr = first;
  ++next;

  if (next == curr)
  {
    // The vertex has a single incident edge. Check if its associated
    // curve equals seg next to vp.
    if (seg_dir_right && curr->direction() == ARR_RIGHT_TO_LEFT)
    {
      // Both curves are defined to the right of vp:
      equal_curr =
        (m_traits->compare_y_at_x_right_2_object() (curr->curve(), seg, vp) ==
         EQUAL);
    }
    else if (! seg_dir_right && curr->direction() == ARR_LEFT_TO_RIGHT)
    {
      // Both curves are defined to the left of vp:
      equal_curr =
        (m_traits->compare_y_at_x_left_2_object() (curr->curve(), seg, vp) ==
         EQUAL);
    }

    // In case the curves are not equal, just return the incident face of
    // the single halfegde (note that this is also the incident face of its
    // twin, as v is the tip of an "antenna").
    if (! equal_curr)
    {
      CGAL_assertion (curr->face() == curr->twin()->face());
      return (CGAL::make_object (curr->face()));
    }
  }
  else
  {
    // Traverse the halfedges around v until we find the pair of adjacent
    // halfedges such as seg is located clockwise in between them.
    typename Traits_adaptor_2::Is_between_cw_2  is_between_cw =
      m_traits->is_between_cw_2_object();
    bool                                        eq_curr, eq_next;

    while (! is_between_cw (seg, seg_dir_right, curr->curve(),
                            (curr->direction() == ARR_RIGHT_TO_LEFT),
                            next->curve(), (next->direction() == ARR_RIGHT_TO_LEFT),
                            vp, eq_curr, eq_next))
    {
      // Break the loop if seg equals one of the halfegdes next to v.
      if (eq_curr)
      {
        equal_curr = true;
        break;
      }

      if (eq_next)
      {
        curr = next;
        equal_curr = true;
        break;
      }

      // Move to the next pair of incident halfedges.
      curr = next;
      ++next;

      // Guard for an infinitive loop, in case we have completed a full
      // traversal around v without locating a place for seg.
      if (curr == first)
      {
        CGAL_assertion_msg
          (false,
           "Completed a full cycle around v without locating seg.");
        return Object();
      }
    }

    // In case seg is not equal to curr's curve, just return the incident face
    // of the halfegde we have located.
    if (! equal_curr)
      return (CGAL::make_object (curr->face()));
  }

  // If we reached here, seg overlaps the curve associated with curr next to
  // the vertex v. We first check if p equals the other end-vertex of this
  // halfedge.
  if (m_traits->equal_2_object() (p, curr->source()->point())) 
  {
    // In this case p equals the source point of the edge.
    return (CGAL::make_object (curr->source()));
  }

  // Check whether p lies on the curve associated with the edge.
  if (m_traits->is_in_x_range_2_object() (curr->curve(), p) && 
      m_traits->compare_y_at_x_2_object() (p, curr->curve()) == EQUAL)
  {
    // p is located on the interior of the edge.
    Halfedge_const_handle   he = curr;
    return (CGAL::make_object (he));
  }

  // In this case, the source vertex of the current edge is closer
  // to the query point p.
  new_vertex = true;
  return (CGAL::make_object (curr->source()));
}

//-----------------------------------------------------------------------------
// Walk from the edge to the query point.
//
template <class Arr, class Gen>
Object Arr_landmarks_point_location<Arr, Gen>::_walk_from_edge
    (Halfedge_const_handle eh,
     const Point_2 & np,
     const Point_2 & p,
     Halfedge_set& crossed_edges) const
{
  CGAL_assertion_msg (! eh->is_fictitious(),
                      "_walk_from_edge() from a fictitious edge.");

  const X_monotone_curve_2& cv = eh->curve() ;
  Comparison_result         res;

  // If p equals one of the edge's endpoints, return the vertex
  // that represents this endpoint.
  if (! eh->source()->is_at_infinity() &&
      m_traits->equal_2_object() (p, eh->source()->point()))
  {
    Vertex_const_handle vh = eh->source();
    return (CGAL::make_object(vh));
  }
  if (! eh->target()->is_at_infinity() &&
      m_traits->equal_2_object() (p, eh->target()->point()))
  {
    Vertex_const_handle vh = eh->target();
    return (CGAL::make_object(vh));
  }

  // Create an initial set of edges that have been crossed, which currently
  // contains only the halfedge we are currently on (and its twin).
  crossed_edges.insert (eh);
  crossed_edges.insert (eh->twin());

  // Check whether p is in the x-range of the edge.
  if (m_traits->is_in_x_range_2_object()(cv, p))
  {
    // If p is in eh's x_range, then we need to check if it is above or below
    // it, so we can orient the halfedge eh accordingly, such that it will be
    // incident to the face that is most likely to contain p.
    res = m_traits->compare_y_at_x_2_object()(p, cv);

    switch (res) { 
      case EQUAL:
        // The edge contains p in its interior:
        return (CGAL::make_object (eh));
      case LARGER:
        // p is above cv: the edge must be oriented from left to right.
        if (eh->direction() == ARR_RIGHT_TO_LEFT)
          eh = eh->twin();
        break;
      case SMALLER:
        // p is below cv: the edge must be oriented from right to left.
        if (eh->direction() == ARR_LEFT_TO_RIGHT)
          eh = eh->twin();
        break;
    }

    // Now walk in the incdient face of eh toward p.
    return (_walk_from_face (eh->face(), np, p, crossed_edges));
  }

  // In this case p is in not in the x-range of eh. We select the end
  // vertex of eh that lies closer to p, and walk from this vertex.
  Vertex_const_handle  vh = eh->source();

  if (vh->is_at_infinity())
  {
    vh = eh->target();
    CGAL_assertion (! vh->is_at_infinity());
  }
  else if (! eh->target()->is_at_infinity())
  {
    res = m_traits->compare_xy_2_object() (p, eh->source()->point());

    if ((eh->direction() == ARR_LEFT_TO_RIGHT && res == LARGER) ||
        (eh->direction() == ARR_RIGHT_TO_LEFT && res == SMALLER))
      vh = eh->target();
  }

  // Walk from the closer vertex toward the query point p.
  return (_walk_from_vertex(vh, p, crossed_edges));
}

//-----------------------------------------------------------------------------
// Walk from the given face to the query point.
//
template <class Arr, class Gen>
Object Arr_landmarks_point_location<Arr, Gen>::_walk_from_face
    (Face_const_handle face,
     const Point_2& np,
     const Point_2& p,
     Halfedge_set& crossed_edges) const
{
  // Construct an x-monotone curve connecting the nearest landmark point np
  // to the query point p and check which CCB intersects this segment.
  X_monotone_curve_2             seg = 
    m_traits->construct_x_monotone_curve_2_object()(np, p);
  const bool                     p_is_left =
    (m_traits->compare_xy_2_object()(np, p) == LARGER);

  Inner_ccb_const_iterator       inner_ccb_iter;
  Outer_ccb_const_iterator       outer_ccb_iter;
  const Halfedge_const_handle    invalid_he;
  Halfedge_const_handle          he;
  Face_const_handle              new_face;
  bool                           is_on_edge;
  bool                           is_target;

  do
  {
    // Check whether p lies inside the current face (including its holes):
    if (p_arr->topology_traits()->is_in_face (&(*face), p, NULL))
    {
      // We know that p is located inside the current face, and we check
      // whether it lies inside one of its holes (or on the boundary of
      // its holes).
      new_face = face;
      for (inner_ccb_iter = face->inner_ccbs_begin();
           inner_ccb_iter != face->inner_ccbs_end(); ++inner_ccb_iter)
      {
        he = _intersection_with_ccb (*inner_ccb_iter,
                                     seg, p, p_is_left,
                                     crossed_edges,
                                     is_on_edge, is_target);

        if (he != invalid_he)
        {
          // Check if the query point is located on a vertex or on an edge.
          if (is_target)
            return (CGAL::make_object (he->target()));
          else if (is_on_edge)
            return (CGAL::make_object (he));

          // Otherwise, cross over he to the incident face of its twin.
          if (face != he->twin()->face())
          {
            new_face = he->twin()->face();
            break;
          }
        }
      }

      // Check if we found a new face (hole) containing p. If not, the current
      // face contains p.
      if (new_face == face)
          return (CGAL::make_object (face));

      // Continue from the new face (hole).
      face = new_face;
    }
    else
    {
      // We know that p is not located inside the current face. We therefore
      // look for an edge on its outer boundary that intersects seg.
      new_face = face;
      for (outer_ccb_iter = face->outer_ccbs_begin();
           outer_ccb_iter != face->outer_ccbs_end(); ++outer_ccb_iter)
      {
        he = _intersection_with_ccb (*outer_ccb_iter,
                                     seg, p, p_is_left,
                                     crossed_edges,
                                     is_on_edge, is_target);

        if (he != invalid_he)
        {
          // Check if the query point is located on a vertex or on an edge.
          if (is_target)
            return (CGAL::make_object (he->target()));
          else if (is_on_edge)
            return (CGAL::make_object (he));

          // Otherwise, cross over he to the incident face of its twin.
          if (face != he->twin()->face())
          {
            new_face = he->twin()->face();
            break;
          }
        }
      }

      // Continue from the new face.
      CGAL_assertion (new_face != face);
      face = new_face;
    }

  } while (true);

  // We should never reach here:
  CGAL_error();
  return Object();
}

//-----------------------------------------------------------------------------
// Find a halfedge on the given CCB that intersects the given x-monotone
// curve, connecting the current landmark to the query point.
//
template <class Arr, class Gen>
typename Arr_landmarks_point_location<Arr, Gen>::Halfedge_const_handle
Arr_landmarks_point_location<Arr, Gen>::_intersection_with_ccb
    (Ccb_halfedge_const_circulator circ,
     const X_monotone_curve_2& seg,
     const Point_2& p, bool p_is_left,
     Halfedge_set& crossed_edges,
     bool& is_on_edge,
     bool& is_target) const
{
  is_on_edge = false;
  is_target = false;

  // Go over the CCB.
  typename Traits_adaptor_2::Is_in_x_range_2    is_in_x_range = 
    m_traits->is_in_x_range_2_object();
  Ccb_halfedge_const_circulator                 curr = circ;
  const Halfedge_const_handle                   invalid_he;
  Halfedge_const_handle                         he;

  do
  {
    he = curr;

    // Skip fictitious halfedges.
    if (he->is_fictitious())
    {
      ++curr;
      continue;
    }

    // Check if we have already crossed the current halfedge (or its twin).
    // If so, we do not cross it again.
    if (crossed_edges.count (he) != 0)
    {
      ++curr;
      continue;
    }

    // Check if the x-range of the curve associated with the current edge
    // does not overlap the x-range of seg, the two curves cannot intersect.
    if (! is_in_x_range (he->curve(), seg))
    {
      ++curr;
      continue;
    }

    // Check whether the current curve intersects seg an odd number of times.
    if (_have_odd_intersections (he->curve(), seg, p_is_left,
                                 is_on_edge))
    {
      // Check if the query point lies on the current edge, or whether
      // it lies in its interior.
      if (is_on_edge)
      {
        // Check if p equals one of the edge end-vertices.
        if (! he->target()->is_at_infinity() &&
            m_traits->compare_xy_2_object() (he->target()->point(), p) == EQUAL)
        {
          // p is the target of the current halfedge.
          is_target = true;
        }
        else if (! he->source()->is_at_infinity() &&
                 m_traits->compare_xy_2_object() (he->source()->point(), p) ==
                 EQUAL)
        {
          // Take the twin halfedge, so p equals its target.
          he = he->twin();
          is_target = true;
        }

        // Return the halfedge containing p.
        return (he);
      }

      // Return the halfedge we found, and mark that we have already crossed
      // it (as well as its twin).
      crossed_edges.insert (he);
      crossed_edges.insert (he->twin());

      return (he);
    }

    // Proceed to the next halfedge along the CCB.
    ++curr;
  } while (curr != circ);

  // If we reached here, we did not find any edge intersecting seg.
  return (invalid_he);
}

//-----------------------------------------------------------------------------
// Check whether the given curve intersects a simple segment, which connects
// the current landmark to the query point, an odd number of times.
//
template <class Arr, class Gen>
bool Arr_landmarks_point_location<Arr, Gen>::_have_odd_intersections
    (const X_monotone_curve_2& cv,
     const X_monotone_curve_2& seg,
     bool p_is_left,
     bool& p_on_curve) const
{
  p_on_curve = false;

  // Use the left and right endpoints of the segment.
  const Point_2&  seg_left = m_traits->construct_min_vertex_2_object() (seg);
  const Point_2&  seg_right = m_traits->construct_max_vertex_2_object() (seg);

  // Check the boundary conditions of the left and right curve ends of cv.
  const Boundary_type  bx_l = m_traits->boundary_in_x_2_object() (cv, ARR_MIN_END);
  const Boundary_type  by_l = m_traits->boundary_in_y_2_object() (cv, ARR_MIN_END);
  const Boundary_type  bx_r = m_traits->boundary_in_x_2_object() (cv, ARR_MAX_END);
  const Boundary_type  by_r = m_traits->boundary_in_y_2_object() (cv, ARR_MAX_END);

  // Check if the overlapping x-range of the two curves is trivial.
  // In this case, they cannot cross.
  if (bx_l != MINUS_INFINITY && bx_l != PLUS_INFINITY &&
      by_l != MINUS_INFINITY && by_l != PLUS_INFINITY)
  {
    // Check if the left endpoint of cv has the same x-coordinate as the
    // right endpoint of seg.
    if (m_traits->compare_x_2_object()
        (m_traits->construct_min_vertex_2_object() (cv), seg_right) == EQUAL)
    {
      if (! p_is_left &&
          m_traits->compare_xy_2_object()
          (m_traits->construct_min_vertex_2_object() (cv), seg_right) == EQUAL)
      {
        p_on_curve = true;
        return (true);
      }
      else if (m_traits->is_vertical_2_object()(seg))
      {
        // Special treatment for vertical segments.
        Comparison_result   res_l =
          m_traits->compare_y_at_x_2_object() (seg_left, cv);
        Comparison_result   res_r =
          m_traits->compare_y_at_x_2_object() (seg_right, cv);

        if ((p_is_left && res_l == EQUAL) ||
            (! p_is_left && res_r == EQUAL))
        {
          p_on_curve = true;
          return (true);
        }

        return (res_l != res_r);
      }

      return (false);
    }
  }
  
  if (bx_r != MINUS_INFINITY && bx_r != PLUS_INFINITY &&
      by_r != MINUS_INFINITY && by_r != PLUS_INFINITY)
  {
    // Check if the right endpoint of cv has the same x-coordinate as the
    // left endpoint of seg.
    if (m_traits->compare_x_2_object()
        (m_traits->construct_max_vertex_2_object() (cv), seg_left) == EQUAL)
    {
      if (p_is_left &&
          m_traits->compare_xy_2_object()
          (m_traits->construct_max_vertex_2_object() (cv), seg_left) == EQUAL)
      {
        p_on_curve = true;
        return (true);
      }
      else if (m_traits->is_vertical_2_object()(seg))
      {
        // Special treatment for vertical segments.
        Comparison_result   res_l =
          m_traits->compare_y_at_x_2_object() (seg_left, cv);
        Comparison_result   res_r =
          m_traits->compare_y_at_x_2_object() (seg_right, cv);

        if ((p_is_left && res_l == EQUAL) ||
            (! p_is_left && res_r == EQUAL))
        {
          p_on_curve = true;
          return (true);
        }

        return (res_l != res_r);
      }

      return (false);
    }
  }

  // Compare the two left ends of cv and seg.
  Comparison_result    left_res;

  if (CGAL::sign (bx_l) == CGAL::NEGATIVE)
  {
    // The left end of cv lies to the left of seg_left:
    // Compare this point to cv.
    left_res = m_traits->compare_y_at_x_2_object() (seg_left, cv);
  }
  else if (CGAL::sign (bx_l) == CGAL::POSITIVE)
  {
    CGAL_assertion (bx_l != PLUS_INFINITY);

    // The left end of cv lies to the right of seg_left.
    // Compare the left endpoint of cv to seg.
    left_res = m_traits->compare_y_at_x_2_object()
        (m_traits->construct_min_vertex_2_object() (cv), seg);
    left_res = CGAL::opposite (left_res);
  }
  else
  {
    if (by_l == MINUS_INFINITY)
      // The left end of cv is at y = -oo, so cv obviously lies above it.
      left_res = LARGER;
    else if (by_l == PLUS_INFINITY)
      // The left end of cv is at y = +oo, so cv obviously lies below it.
      left_res = SMALLER;
    else
    {
      // In this case cv has a valid left endpoint: Find the rightmost of
      // these two points and compare it to the other curve.
      const Point_2&    cv_left = m_traits->construct_min_vertex_2_object()(cv);
      Comparison_result res =
        m_traits->compare_xy_2_object() (cv_left, seg_left);

      if (res != LARGER)
      {
        left_res = m_traits->compare_y_at_x_2_object() (seg_left, cv);
      }
      else
      {
        left_res = m_traits->compare_y_at_x_2_object() (cv_left, seg);
        left_res = CGAL::opposite (left_res);
      }
    }
  }

  if (left_res == EQUAL)
  {
    if (p_is_left)
    {
      // In this case the query point p, which is the left endpoint of seg,
      // lies on cv.
      p_on_curve = true;
      return (true);
    }

    // Compare the two curves to the right of their common left endpoint.
    left_res = m_traits->compare_y_at_x_right_2_object() (seg, cv, seg_left);

    if (left_res == EQUAL)
    {
        // RWRW: In this case we have an overlap ...
    }
  }

  // Compare the two right ends of cv and seg.
  Comparison_result    right_res;

  if (CGAL::sign (bx_r) == CGAL::POSITIVE)
  {
    // The right end of cv lies to the right of seg_right:
    // Compare this point to cv.
    right_res = m_traits->compare_y_at_x_2_object() (seg_right, cv);
  }
  else if (CGAL::sign (bx_r) == CGAL::NEGATIVE)
  {
    CGAL_assertion (bx_r != MINUS_INFINITY);

    // The right end of cv lies to the left of seg_right.
    // Compare the right endpoint of cv to seg.
    right_res = m_traits->compare_y_at_x_2_object()
        (m_traits->construct_max_vertex_2_object() (cv), seg);
    right_res = CGAL::opposite (right_res);
  }
  else
  {
    if (by_r == MINUS_INFINITY)
      // The right end of cv is at y = -oo, so cv obviously lies above it.
      right_res = LARGER;
    else if (by_r == PLUS_INFINITY)
      // The right end of cv is at y = +oo, so cv obviously lies below it.
      right_res = SMALLER;
    else
    {
      // In this case cv has a valid right endpoint: Find the leftmost of
      // these two points and compare it to the other curve.
      const Point_2& cv_right = m_traits->construct_max_vertex_2_object()(cv);
      Comparison_result res =
        m_traits->compare_xy_2_object() (cv_right, seg_right);

      if (res != SMALLER)
      {
        right_res = m_traits->compare_y_at_x_2_object() (seg_right, cv);
      }
      else
      {
        right_res = m_traits->compare_y_at_x_2_object() (cv_right, seg);
        right_res = CGAL::opposite (right_res);
      }
    }
  }

  if (right_res == EQUAL)
  {
    if (! p_is_left)
    {
      // In this case the query point p, which is the right endpoint of seg,
      // lies on cv.
      p_on_curve = true;
      return (true);
    }

    // Compare the two curves to the left of their common right endpoint.
    right_res = m_traits->compare_y_at_x_left_2_object() (seg, cv, seg_right);

    if (right_res == EQUAL)
    {
        // RWRW: In this case we have an overlap ...
    }
  }

  // The two curves intersect an odd number of times if the comparison
  // results at the two ends are not the same (this indicates that they
  // switch positions).
  return (left_res != right_res);
}

CGAL_END_NAMESPACE

#endif
