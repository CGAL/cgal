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
   * given point hits (not inculding isolated vertices).
   * \param p The query point.
   * \param shoot_up Indicates whether the ray is directed upward or downward.
   * \return An object representing the arrangement feature the ray hits.
   *         This object is either a Halfedge_const_handle,
   *         a Vertex_const_handle or an empty object.
   */
  Object _base_vertical_ray_shoot (const Point_2& p, bool shoot_up) const;

  /*!
   * Locate the arrangement feature which a vertical ray emanating from the
   * given point hits, considering isolated vertices.
   * \param p The query point.
   * \param shoot_up Indicates whether the ray is directed upward or downward.
   * \return An object representing the arrangement feature the ray hits.
   *         This object is either a Halfedge_const_handle,
   *         a Vertex_const_handle or an empty object.
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

CGAL_END_NAMESPACE

// The member-function definitions can be found under:
#include <CGAL/Arr_point_location/Arr_naive_point_location_functions.h>

#endif
