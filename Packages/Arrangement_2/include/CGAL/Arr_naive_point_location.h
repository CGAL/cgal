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
#ifndef CGAL_ARR_NAIVE_POINT_LOCATION_H
#define CGAL_ARR_NAIVE_POINT_LOCATION_H

/*! \file
 * Definition of the Arr_naive_point_location<Arrangement> template.
 */

#include <CGAL/Arrangement_2/Arr_traits_wrapper_2.h>

CGAL_BEGIN_NAMESPACE

/*! \class
 * A class that answers point-location and vertical ray-shooting queries
 * on a planar arrangement using a naive algorithm.
 * The Arrangement parameter corresponds to an arrangement instantiation.
 */
template <class Arrangement_>
class Arr_naive_point_location
{
public:

  typedef Arrangement_                         Arrangement_2;
  typedef typename Arrangement_2::Traits_2     Traits_2;

  typedef typename Arrangement_2::Vertex_const_handle   Vertex_const_handle;
  typedef typename Arrangement_2::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Arrangement_2::Face_const_handle     Face_const_handle;

  typedef typename Traits_2::Point_2            Point_2;
  typedef typename Traits_2::X_monotone_curve_2 X_monotone_curve_2;

protected:

  typedef Arr_traits_basic_wrapper_2<Traits_2>  Traits_wrapper_2;

  // Data members:
  const Arrangement_2     *p_arr;     // The associated arrangement.
  const Traits_wrapper_2  *traits;    // Its associated traits object.
        
public:

  /*! Default constructor. */
  Arr_naive_point_location () : 
    p_arr (NULL),
    traits (NULL)
  {}
        
  /*! Constructor given an arrangement. */
  Arr_naive_point_location (const Arrangement_2& arr) :
    p_arr (&arr)
  {
    traits = static_cast<const Traits_wrapper_2*> (p_arr->get_traits());
  }
        
  /*! Attach an arrangement object. */
  void init (const Arrangement_2& arr) 
  {
    p_arr = &arr;
    traits = static_cast<const Traits_wrapper_2*> (p_arr->get_traits());
  }
        
  /*!
   * Locate the arrangement feature containing the given point.
   * \param p The query point.
   * \return An object representing the arrangement feature containing the
   *         query point. This object is either a Face_const_handle or a
   *         Halfedge_const_handle or a Vertex_const_handle.
   */
  Object locate (const Point_2& p) const;

  /*!
   * Locate the arrangement feature which a upward vertical ray emanating from
   * the given point hits.
   * \param p The query point.
   * \return An object representing the arrangement feature the ray hits.
   *         This object is either an empty object or a
   *         Halfedge_const_handle or a Vertex_const_handle.
   */
  Object ray_shoot_up (const Point_2& p) const
  {
    return (_vertical_ray_shoot (p, true));
  }

  /*!
   * Locate the arrangement feature which a downward vertical ray emanating
   * from the given point hits.
   * \param p The query point.
   * \return An object representing the arrangement feature the ray hits.
   *         This object is either an empty object or a
   *         Halfedge_const_handle or a Vertex_const_handle.
   */
  Object ray_shoot_down (const Point_2& p) const
  {
    return (_vertical_ray_shoot (p, false));
  }

protected:

  /*!
   * Locate the arrangement feature which a vertical ray emanating from the
   * given point hits.
   * \param p The query point.
   * \param shoot_up Indicates whether the ray is directed upward or downward.
   * \return An object representing the arrangement feature the ray hits.
   *         This object is either a Face_const_handle or a
   *         Halfedge_const_handle or a Vertex_const_handle.
   */
  Object _vertical_ray_shoot (const Point_2& p, bool shoot_up) const;

  /*!
   * Find the first halfedge with a given source vertex, when going clockwise
   * from "6 o'clock" around this vertex.
   * \param v The given vertex.
   * \return The first halfedge.
   */
  Halfedge_const_handle _first_around_vertex (Vertex_const_handle v) const;
};

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

  for (vit = p_arr->vertices_begin(); vit != p_arr->vertices_end(); ++vit)
  {
    if (equal (p, (*vit).point()))
      return (make_object (*vit));
  }

  // Go over arrangement halfedges and check whether one of them contains
  // the query point in its interior.
  typename Traits_wrapper_2::Is_in_x_range_2    is_in_x_range = 
                                            traits->is_in_x_range_2_object();
  typename Traits_wrapper_2::Compare_y_at_x_2   compare_y_at_x = 
                                            traits->compare_y_at_x_2_object();
  typename Arrangement::Edge_const_iterator     eit;

  for (eit = p_arr->edges_begin(); eit != p_arr->edges_end(); ++eit)
  {
    if (is_in_x_range ((*eit).curve(), p) &&
        compare_y_at_x (p, (*eit).curve()) == EQUAL)
    {
      return (make_object (*eit));
    }
  }

  // Shoot a vertical ray from the query point.
  Object   obj = _vertical_ray_shoot (p, true);

  if (obj.is_empty())
  {
    // Return the unbounded face.
    return (make_object (p_arr->unbounded_face()));
  }

  // The ray shooting returned either a vertex of a halfedge.
  typename Arrangement::Halfedge_const_handle   h;

  if (assign (h, obj))
  {
    // Make sure that the edge is directed from right to left, so that p
    // (which lies below it) is contained in its incident face. If necessary,
    // we take the twin halfedge.
    if (traits->compare_xy_2_object() (h.source().point(), 
                                       h.target().point()) == SMALLER)
    {
      h = h.twin();
    }

    // Return the incident face.
    return (make_object (h.face()));
  }

  // In case the ray-shooting returned a vertex, we have to locate the first
  // halfedge whose source vertex is v, rotating clockwise around the vertex
  // from "6 o'clock", and to return its incident face. 
  typename Arrangement::Vertex_const_handle     v;

  if (! assign (v, obj))
  {
    CGAL_assertion (false);
    return Object();
  }

  h = _first_around_vertex (v);
  return (make_object (h.face()));
}

//-----------------------------------------------------------------------------
// Locate the arrangement feature which a vertical ray emanating from the
// given point hits.
//
template <class Arrangement>
Object Arr_naive_point_location<Arrangement>::_vertical_ray_shoot
    (const Point_2& p,
     bool shoot_up) const
{
  // set the flags for comparison acording to the ray direction.
  Comparison_result res;
  Comparison_result point_above_under;
  Comparison_result curve_above_under;
  bool              in_x_range;

  if (shoot_up)
  {
    point_above_under = SMALLER;
    curve_above_under = LARGER;
  }
  else
  {
    point_above_under = LARGER;
    curve_above_under = SMALLER;
  }

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
  bool                                         found = false;
  
  while (eit != p_arr->edges_end())
  {
    // Determine whether p is in the x-range of the curve and above or below it
    // (according to the direction of the shoot).
    in_x_range = is_in_x_range ((*eit).curve(), p);
    if (in_x_range)
      res = compare_y_at_x (p, (*eit).curve());

    if (in_x_range && res == point_above_under) 
    {
      if (!found)
      {
        // If no other x-monotone curve containing p in its x-range has been
        // found yet, take the current one as the vertically closest to p.
        closest_edge = *eit;
        found = true;
      }
      else
      {
        // Compare with the vertically closest cure so far and detemine the
        // curve closest to p. Note that the two curves do not intersect
        // in their interiors.
        if (compare_y_position (closest_edge.curve(),
                                (*eit).curve()) == curve_above_under)
        {
          closest_edge = *eit;
        }
      }
    }

    if (in_x_range && res == EQUAL &&
        is_vertical((*eit).curve()))
    {
      // The vertical ray overlaps an existing vertical edge containing p.
      // In this case simply return this edge.
      return (make_object (*eit));
    }

    // Move to the next edge.
    ++eit;
  }
  
  // If we have not found any edge above p, we return an empty object.
  if (!found)
    return Object();

  // If one of the closest edge's end vertices has the same x-coordinate
  // as the query point, return this vertex.
  if (traits->compare_x_2_object() (closest_edge.source().point(), 
                                    p) == EQUAL)
  {
    return (make_object (closest_edge.source()));
  }
  else if (traits->compare_x_2_object() (closest_edge.target().point(), 
                                         p) == EQUAL)
  {
    return (make_object (closest_edge.target()));
  }

  // Otherwise, return the closest edge.
  return (make_object (closest_edge));
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
    v.incident_halfedges();
  typename Arrangement::Halfedge_around_vertex_const_circulator curr = first;

  do 
  {
    // Check whether the current halfedge is defined to the left or to the
    // right of the given vertex.
    if (compare_xy ((*curr).source().point(), v.point()) == SMALLER)
    {
      // The curve associated with the current halfedge is defined to the left
      // of v.
      if (lowest_left == invalid_handle ||
          compare_y_at_x_left ((*curr).curve(),
                               lowest_left.curve(), 
                               v.point()) == SMALLER)
      {
        lowest_left = *curr;
      }
    }
    else
    {
      // The curve associated with the current halfedge is defined to the right
      // of v.
      if (top_right == invalid_handle ||
          compare_y_at_x_right ((*curr).curve(),
                                top_right.curve(), 
                                v.point()) == LARGER)
      {
        top_right = *curr;
      }
    }

    curr++;
  } while (curr != first);

  // The first halfedge we encounter is the lowest to the left, but if there
  // is no edge to the left, we first encounter the topmost halfedge to the 
  // right. Note that as the halfedge we located has v as its target, we now
  // have to return its twin.
  if (lowest_left != invalid_handle)
    return (lowest_left.twin());
  else
    return (top_right.twin());
}

CGAL_END_NAMESPACE

#endif //CGAL_PM_NAIVE_POINT_LOCATION_H

