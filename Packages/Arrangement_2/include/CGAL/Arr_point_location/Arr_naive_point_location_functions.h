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
// Author(s)     : Ron Wein   <wein@post.tau.ac.il>
//                 (based on old version by Eyal Flato)
#ifndef CGAL_ARR_NAIVE_POINT_LOCATION_FUNCTIONS_H
#define CGAL_ARR_NAIVE_POINT_LOCATION_FUNCTIONS_H

/*! \file
 * Member-function definitions for the Arr_naive_point_location<Arrangement>
 * class.
 */

CGAL_BEGIN_NAMESPACE

//-----------------------------------------------------------------------------
// Locate the arrangement feature containing the given point.
//
template <class Arrangement>
Object Arr_naive_point_location<Arrangement>::locate (const Point_2& p) const
{
  // Go over the arrangement vertices and check whether one of them equals
  // the query point.
  typename Traits_wrapper_2::Equal_2            equal = 
                                            traits->equal_2_object();
  typename Arrangement::Vertex_const_iterator   vit;
  typename Arrangement::Vertex_const_handle     vh;

  for (vit = p_arr->vertices_begin(); vit != p_arr->vertices_end(); ++vit)
  {
    vh = vit;
    if (equal (p, vh->point()))
      return (CGAL::make_object (vh));
  }

  // Go over arrangement halfedges and check whether one of them contains
  // the query point in its interior.
  typename Traits_wrapper_2::Is_in_x_range_2    is_in_x_range = 
                                            traits->is_in_x_range_2_object();
  typename Traits_wrapper_2::Compare_y_at_x_2   compare_y_at_x = 
                                            traits->compare_y_at_x_2_object();
  typename Arrangement::Edge_const_iterator     eit;
  typename Arrangement::Halfedge_const_handle   hh;

  for (eit = p_arr->edges_begin(); eit != p_arr->edges_end(); ++eit)
  {
    hh = eit;
    if (is_in_x_range (hh->curve(), p) &&
        compare_y_at_x (p, hh->curve()) == EQUAL)
    {
      return (CGAL::make_object (hh));
    }
  }

  // Shoot a vertical ray from the query point.
  Object   obj = _base_vertical_ray_shoot (p, true);

  if (obj.is_empty())
  {
    // Return the unbounded face.
    return (CGAL::make_object (p_arr->unbounded_face()));
  }

  // The ray shooting returned either a vertex of a halfedge.
  const typename Arrangement::Vertex_const_handle    *p_vh;
  const typename Arrangement::Halfedge_const_handle  *p_hh;

  p_hh = object_cast<typename Arrangement::Halfedge_const_handle> (&obj);
  if (p_hh != NULL)
  {
    // Make sure that the edge is directed from right to left, so that p
    // (which lies below it) is contained in its incident face. If necessary,
    // we take the twin halfedge.
    if (traits->compare_xy_2_object() ((*p_hh)->source()->point(), 
                                       (*p_hh)->target()->point()) == SMALLER)
    {
      hh = (*p_hh)->twin();
    }
    else
    {
      hh = *p_hh;
    }

    // Return the incident face.
    return (CGAL::make_object (hh->face()));
  }

  // In case the ray-shooting returned a vertex, we have to locate the first
  // halfedge whose source vertex is v, rotating clockwise around the vertex
  // from "6 o'clock", and to return its incident face. 
  p_vh = object_cast<typename Arrangement::Vertex_const_handle> (&obj);
  CGAL_assertion (p_vh != NULL);

  hh = _first_around_vertex (*p_vh);
  return (CGAL::make_object (hh->face()));
}

//-----------------------------------------------------------------------------
// Locate the arrangement feature which a vertical ray emanating from the
// given point hits (not inculding isolated vertices).
//
template <class Arrangement>
Object Arr_naive_point_location<Arrangement>::_base_vertical_ray_shoot
    (const Point_2& p,
     bool shoot_up) const
{
  // Set the results for comparison according to the ray direction.
  const Comparison_result point_above_under = (shoot_up ? SMALLER : LARGER);
  const Comparison_result curve_above_under = (shoot_up ? LARGER : SMALLER);

  // Go over all halfedges in the arrangement.
  typename Traits_wrapper_2::Is_in_x_range_2      is_in_x_range =
                                        traits->is_in_x_range_2_object();
  typename Traits_wrapper_2::Compare_y_at_x_2     compare_y_at_x =
                                        traits->compare_y_at_x_2_object();
  typename Traits_wrapper_2::Is_vertical_2        is_vertical =
                                        traits->is_vertical_2_object();
  typename Traits_wrapper_2::Compare_y_position_2 compare_y_position =
                                        traits->compare_y_position_2_object();

  typename Arrangement::Edge_const_iterator    eit = p_arr->edges_begin();
  typename Arrangement::Halfedge_const_handle  closest_edge;
  Comparison_result                            res;
  bool                                         in_x_range;
  bool                                         found = false;

  while (eit != p_arr->edges_end())
  {
    // Determine whether p is in the x-range of the curve and above or below it
    // (according to the direction of the shoot).
    in_x_range = is_in_x_range (eit->curve(), p);
    if (in_x_range)
      res = compare_y_at_x (p, eit->curve());

    if (in_x_range && res == point_above_under)
    {
      if (!found)
      {
        // If no other x-monotone curve containing p in its x-range has been
        // found yet, take the current one as the vertically closest to p.
        closest_edge = eit;
        found = true;
      }
      else
      {
        // Compare with the vertically closest cure so far and detemine the
        // curve closest to p. Note that the two curves do not intersect
        // in their interiors.
        if (compare_y_position (closest_edge->curve(),
                                eit->curve()) == curve_above_under)
        {
          closest_edge = eit;
        }
      }
    }

    if (in_x_range && res == EQUAL &&
        is_vertical(eit->curve()))
    {
      // The vertical ray overlaps an existing vertical edge containing p.
      // In this case simply return this edge.
      return (CGAL::make_object (eit));
    }

    // Move to the next edge.
    ++eit;
  }

  // If we have not found any edge above p, we return an empty object.
  if (!found)
    return Object();

  // If one of the closest edge's end vertices has the same x-coordinate
  // as the query point, return this vertex.
  if (traits->compare_x_2_object() (closest_edge->source()->point(),
                                    p) == EQUAL)
  {
    return (CGAL::make_object (closest_edge->source()));
  }
  else if (traits->compare_x_2_object() (closest_edge->target()->point(),
                                         p) == EQUAL)
  {
    return (CGAL::make_object (closest_edge->target()));
  }

  // Otherwise, return the closest edge.
  return (CGAL::make_object (closest_edge));
}

//-----------------------------------------------------------------------------
// Locate the arrangement feature which a vertical ray emanating from the
// given point hits, considering isolated vertices.
//
template <class Arrangement>
Object Arr_naive_point_location<Arrangement>::_vertical_ray_shoot
    (const Point_2& p,
     bool shoot_up) const
{
  // Locate the arrangement feature which a vertical ray emanating from the
  // given point hits, when not considering the isolated vertices.
  // This feature may not exist, or be either a vertex of a halfedge.
  Object                 obj = _base_vertical_ray_shoot (p, shoot_up);
  Vertex_const_handle    closest_v;
  Halfedge_const_handle  closest_he;

  enum {NAIVE_PL_NONE, NAIVE_PL_VERTEX, NAIVE_PL_HALFEDGE}  type;

  if (obj.is_empty())
  {
    type = NAIVE_PL_NONE;
  }
  else
  {
    const Vertex_const_handle *p_vh = object_cast<Vertex_const_handle> (&obj);
    
    if (p_vh != NULL)
    {
      closest_v = *p_vh;
      type = NAIVE_PL_VERTEX;
    }
    else
    {
      closest_he = object_cast<Halfedge_const_handle> (obj);
      type = NAIVE_PL_HALFEDGE;
    }
  }

  // Set the result for comparison according to the ray direction.
  const Comparison_result point_above_under = (shoot_up ? SMALLER : LARGER);

  // Go over all isolated vertices in the arrangement.
  typename Traits_wrapper_2::Compare_x_2          compare_x =
                                        traits->compare_x_2_object();
  typename Traits_wrapper_2::Compare_xy_2         compare_xy =
                                        traits->compare_xy_2_object();
  typename Traits_wrapper_2::Compare_y_at_x_2     compare_y_at_x =
                                        traits->compare_y_at_x_2_object();

  typename Arrangement::Vertex_const_iterator  vit;
  Vertex_const_handle                          vh;

  for (vit = p_arr->vertices_begin(); vit != p_arr->vertices_end(); ++vit)
  {
    vh = vit;
    if (! vh->is_isolated())
      continue;

    // The current isolated vertex should have the same x-coordinate as the
    // query point in order to be below or above it.
    if (compare_x (p, vh->point()) != EQUAL)
      continue;

    // Make sure the isolated vertex is above the query point (if we shoot up)
    // or below it (if we shoot down).
    if (compare_xy (p, vh->point()) != point_above_under)
      continue;

    // Check if the isolated vertex is closer to p than the current closest
    // object.
    if ((type == NAIVE_PL_NONE) ||
        (type == NAIVE_PL_VERTEX &&
         compare_xy (vh->point(), closest_v->point()) == point_above_under) ||
        (type == NAIVE_PL_HALFEDGE &&
         compare_y_at_x (vh->point(), closest_he->curve()) == 
                                                         point_above_under))
    {
      closest_v = vh;
      type = NAIVE_PL_VERTEX;
    }
  }

  // Set back the result according to its type.
  if (type == NAIVE_PL_NONE)
    return Object();
  else if (type == NAIVE_PL_VERTEX)
    return (CGAL::make_object (closest_v));
  else
    return (CGAL::make_object (closest_he));
}

//-----------------------------------------------------------------------------
// Find the first halfedge with a given source vertex, when going clockwise
// from "6 o'clock" around this vertex.
//
template <class Arrangement>
typename Arr_naive_point_location<Arrangement>::Halfedge_const_handle
Arr_naive_point_location<Arrangement>::_first_around_vertex
    (Vertex_const_handle v) const
{
  // Travrse the incident halfedges of the current vertex and locate the
  // lowest one to its left and the topmost to its right.
  typename Traits_wrapper_2::Compare_xy_2           compare_xy = 
                                      traits->compare_xy_2_object();
  typename Traits_wrapper_2::Compare_y_at_x_right_2 compare_y_at_x_right =
                                      traits->compare_y_at_x_right_2_object();
  typename Traits_wrapper_2::Compare_y_at_x_left_2  compare_y_at_x_left =
                                      traits->compare_y_at_x_left_2_object();

  Halfedge_const_handle   invalid_handle;
  Halfedge_const_handle   lowest_left;
  Halfedge_const_handle   top_right;

  typename Arrangement::Halfedge_around_vertex_const_circulator first = 
    v->incident_halfedges();
  typename Arrangement::Halfedge_around_vertex_const_circulator curr = first;

  do 
  {
    // Check whether the current halfedge is defined to the left or to the
    // right of the given vertex.
    if (compare_xy (curr->source()->point(), v->point()) == SMALLER)
    {
      // The curve associated with the current halfedge is defined to the left
      // of v.
      if (lowest_left == invalid_handle ||
          compare_y_at_x_left (curr->curve(),
                               lowest_left->curve(), 
                               v->point()) == SMALLER)
      {
        lowest_left = curr;
      }
    }
    else
    {
      // The curve associated with the current halfedge is defined to the right
      // of v.
      if (top_right == invalid_handle ||
          compare_y_at_x_right (curr->curve(),
                                top_right->curve(), 
                                v->point()) == LARGER)
      {
        top_right = curr;
      }
    }

    curr++;
  } while (curr != first);

  // The first halfedge we encounter is the lowest to the left, but if there
  // is no edge to the left, we first encounter the topmost halfedge to the 
  // right. Note that as the halfedge we located has v as its target, we now
  // have to return its twin.
  if (lowest_left != invalid_handle)
    return (lowest_left->twin());
  else
    return (top_right->twin());
}

CGAL_END_NAMESPACE

#endif
