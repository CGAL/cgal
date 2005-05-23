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
//                 (based on old version by Oren Nechushtan
//                                      and Iddo Hanniel)
#ifndef CGAL_ARR_WALK_ALONG_LINE_POINT_LOCATION_H
#define CGAL_ARR_WALK_ALONG_LINE_POINT_LOCATION_H

/*! \file
 * Definition of the Arr_walk_along_line_point_location<Arrangement> template.
 */

#include <CGAL/Arrangement_2/Arr_traits_wrapper_2.h>

CGAL_BEGIN_NAMESPACE

/*! \class
 * A class that answers point-location and vertical ray-shooting queries
 * on a planar arrangement by walking on a vertical ray emanating from the
 * query point, going from "infinity" (the unbounded face) until reaching
 * the point.
 * The Arrangement parameter corresponds to an arrangement instantiation.
 */
template <class Arrangement_>
class Arr_walk_along_line_point_location
{
public:

  typedef Arrangement_                          Arrangement_2;
  typedef typename Arrangement_2::Traits_2      Traits_2;

  typedef typename Arrangement_2::Vertex_const_handle   Vertex_const_handle;
  typedef typename Arrangement_2::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Arrangement_2::Face_const_handle     Face_const_handle;

  typedef typename Traits_2::Point_2            Point_2;
  typedef typename Traits_2::X_monotone_curve_2 X_monotone_curve_2;

protected:

  typedef Arr_traits_basic_wrapper_2<Traits_2>  Traits_wrapper_2;
  typedef typename Arrangement_2::Ccb_halfedge_const_circulator
                                                Ccb_halfedge_const_circulator;
  typedef typename Arrangement_2::Holes_const_iterator
                                                Holes_const_iterator;

  // Data members:
  const Arrangement_2     *p_arr;     // The associated arrangement.
  const Traits_wrapper_2  *traits;    // Its associated traits object.
        
public:

  /*! Default constructor. */
  Arr_walk_along_line_point_location () : 
    p_arr (NULL),
    traits (NULL)
  {}
        
  /*! Constructor given an arrangement. */
  Arr_walk_along_line_point_location (const Arrangement_2& arr) :
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
   * \struct An enumeration type specifying how to treat the halfedge returned
   *         by the private functions.
   */
  enum Locate_type
  {
    WPL_VERTEX,            // Take the target vertex of the halfedge.
    WPL_EDGE,              // Take the halfedge itself.
    WPL_FACE,              // Take the face incident to the halfedge.
    WPL_UFACE              // Take the unbounded face.
  };
  
  /*!
   * Walk along a the vertical line enamating from p from the given halfedge.
   * \param p The query point.
   * \param shoot_up Whether we shoot the vertical ray up or down.
   * \param inclusive Indicates whether the vertical ray includes the point p:
   *                  If not (in case of vertical ray shooting) we should find
   *                  the edge or vertex right above (or below) the query point
   *                  p. If it does (in case of point location) and p lies on
   *                  a vertex of on an edge, we return this feature.
   * \param closest_he Input: The closest halfedge to p so far.
   *                   Output: The updated closest halfedge.
   * \param locate_type Input: How to treat the input handle closest_he.
   *                    Output: How to treat the output handle closest_he.
   */
  void _walk_along_line (const Point_2& p,
			 bool shoot_up,
			 bool inclusive,
			 Halfedge_const_handle& he,
			 Locate_type& lt) const;

  /*!
   * Find the closest feature to p (and lying above or below it) along the
   * boundary of the given connected component. 
   * \param p The query point.
   * \param circ A circulator for the halfedges along the connected component
   *             boundary.
   * \param shoot_up Whether we shoot the vertical ray up or down.
   * \param inclusive Indicates whether the vertical ray includes the point p:
   *                  If not (in case of vertical ray shooting) we should find
   *                  the edge or vertex right above (or below) the query point
   *                  p. If it does (in case of point location) and p lies on
   *                  a vertex of on an edge, we return this feature.
   * \param closest_he Input: The closest halfedge to p so far.
   *                   Output: The updated closest halfedge.
   * \param locate_type Input: How to treat the input handle closest_he.
   *                    Output: How to treat the output handle closest_he.
   * \return (true) if the closest_he is valid; (false) otherwise.
   */
  bool _find_closest_feature (const Point_2& p,
			      Ccb_halfedge_const_circulator circ,
			      bool shoot_up, bool inclusive,
			      Halfedge_const_handle& closest_he,
			      Locate_type& locate_type) const;
};

CGAL_END_NAMESPACE

// The member-function definitions can be found under:
#include <CGAL/Arr_point_location/Arr_walk_along_line_pl_functions.h>

#endif
