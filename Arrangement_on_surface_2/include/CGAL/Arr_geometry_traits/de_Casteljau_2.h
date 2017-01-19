// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
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
// Author(s)     : Ron Wein     <wein@post.tau.ac.il>
//                 Iddo Hanniel <iddoh@cs.technion.ac.il>

#ifndef CGAL_DE_CASTELJAU_2_H
#define CGAL_DE_CASTELJAU_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>


/*! \file
 * Template functions for performing generic operations based on the
 * de Casteljau algorithm.
 */

#include <vector>

namespace CGAL {

/*!
 * Bisect the control polygon of a given Bezier curve into the left and right
 * control polygons.
 * \param ctrl_pts_begin The beginning iterator of the control points of the
 *        curve.
 * \param ctrl_pts_end The past-the-end iterator of the control points.
 * \param left_ctrl_pts (out) The control points of the left polygon.
 * \param right_ctrl_pts (out) The control points of the right polygon.
 * Note: Typically you should call this function as follows:
 *       bisect_control_polygon_2(ctrl_pts.begin(), ctrl_pts.end(),
 *                                std::back_inserter(left_pts),
 *                                std::front_inserter(right_pts));
 */
template <class InputIterator,
          class OutputIterLeft, class OutputIterRight>
void
bisect_control_polygon_2(InputIterator ctrl_pts_begin,
                         InputIterator ctrl_pts_end,
                         OutputIterLeft left_ctrl_pts,
                         OutputIterRight right_ctrl_pts)
{
  typedef typename InputIterator::value_type        _Point_2;
  typedef typename Kernel_traits<_Point_2>::Kernel  _Kernel;
  typedef typename _Kernel::Construct_midpoint_2    _Construct_midpoint_2;

  // Grab a local copy of the control points.
  const unsigned int    n_pts = std::distance(ctrl_pts_begin, ctrl_pts_end);
  CGAL_precondition(n_pts != 0);

  std::vector<_Point_2>  vec(n_pts);
  InputIterator          iter;
  unsigned int           i;

  for (iter = ctrl_pts_begin, i = 0; iter != ctrl_pts_end; ++iter, ++i)
    vec[i] = *iter;

  // The first control point goes to the (front of) the right subcurve,
  // while the last control point goes to the (back of) the left subcurve.
  _Kernel                ker;
  _Construct_midpoint_2  midpoint = ker.construct_midpoint_2_object();
  unsigned int           last_index = n_pts - 1;

  *left_ctrl_pts = vec[0];
  ++left_ctrl_pts;

  *right_ctrl_pts = vec[last_index];
  ++right_ctrl_pts;

  while (last_index > 0)
  {
    // Construct (m - 1) control points from the m point we currently have,
    // where the new i'th point is the midpoint between p[i]  and p[i + 1].
    for (i = 0; i < last_index; ++i)
      vec[i] = midpoint(vec[i], vec[i + 1]);

    --last_index;

    // The first control points goes to the (front of) the right subcurve,
    // while the last control point goes to the (back of) the left subcurve.
    *left_ctrl_pts = vec[0];
    ++left_ctrl_pts;

    *right_ctrl_pts = vec[last_index];
    ++right_ctrl_pts;
  }

  return;
}

/*!
 * Evaluate a point on a parametric Bezier curve B(t) using de Casteljau's
 * algorithm.
 * \param ctrl_pts_begin The beginning iterator of the control points of the
 *        curve.
 * \param ctrl_pts_end The past-the-end iterator of the control points.
 * \param t0 The parameter value at the requested point.
 * \return The point B(t0).
 */
template <class InputIterator>
typename InputIterator::value_type point_on_Bezier_curve_2
(InputIterator ctrl_pts_begin,
 InputIterator ctrl_pts_end,
 const typename Kernel_traits<typename InputIterator::value_type>::Kernel::FT&
 t0)
{
  typedef typename InputIterator::value_type        _Point_2;
  typedef typename Kernel_traits<_Point_2>::Kernel  _Kernel;
  typedef typename _Kernel::FT                      _NT;

  // Grab a local copy of the control points.
  const unsigned int     n_pts = std::distance(ctrl_pts_begin, ctrl_pts_end);
  CGAL_precondition(n_pts != 0);

  std::vector<_Point_2>  vec(n_pts);
  InputIterator          iter;
  unsigned int           i;

  for (iter = ctrl_pts_begin, i = 0; iter != ctrl_pts_end; ++iter, ++i)
    vec[i] = *iter;

  // The first control point goes to the (front of) the right subcurve,
  // while the last control point goes to the (back of) the left subcurve.
  const _NT              comp_t0 = _NT(1) - t0;
  unsigned int           last_index = n_pts - 1;

  while (last_index > 0)
  {
    // Construct (m - 1) control points from the m point we currently have,
    // where the new i'th point is given by: (1 - t0)*p[i] + t0*p[i + 1].
    for (i = 0; i < last_index; ++i)
    {
      vec[i] = _Point_2(comp_t0*vec[i].x() + t0*vec[i + 1].x(),
                        comp_t0*vec[i].y() + t0*vec[i + 1].y());
    }

    --last_index;
  }

  // Our auxiliary vector now contains just a single entry, and this entry
  // is the curve value at t0.
  return (vec[0]);
}

/*!
 * Evaluate a point on a parametric Bezier curve B(t) using de Casteljau's
 * algorithm, and also bisect the control polygon of B at this point.
 * \param ctrl_pts_begin The beginning iterator of the control points of the
 *        curve.
 * \param ctrl_pts_end The past-the-end iterator of the control points.
 * \param t0 The parameter value at the requested point.
 * \param left_ctrl_pts (out) The control points of the left polygon.
 * \param right_ctrl_pts (out) The control points of the right polygon.
 * \return The point B(t0).
 */
template <class InputIterator,
          class OutputIterLeft, class OutputIterRight>
typename InputIterator::value_type de_Casteljau_2
(InputIterator ctrl_pts_begin, InputIterator ctrl_pts_end,
 const typename Kernel_traits<typename InputIterator::value_type>::Kernel::FT&
 t0,
 OutputIterLeft left_ctrl_pts, OutputIterRight right_ctrl_pts)
{
  typedef typename InputIterator::value_type        _Point_2;
  typedef typename Kernel_traits<_Point_2>::Kernel  _Kernel;
  typedef typename _Kernel::FT                      _NT;

  // Grab a local copy of the control points.
  const unsigned int     n_pts = std::distance(ctrl_pts_begin, ctrl_pts_end);
  CGAL_precondition(n_pts != 0);

  std::vector<_Point_2>  vec(n_pts);
  InputIterator          iter;
  unsigned int           i;

  for (iter = ctrl_pts_begin, i = 0; iter != ctrl_pts_end; ++iter, ++i)
    vec[i] = *iter;

  // The first control point goes to the (front of) the right subcurve,
  // while the last control point goes to the (back of) the left subcurve.
  const _NT              comp_t0 = _NT(1) - t0; 
  unsigned int           last_index = n_pts - 1;

  *left_ctrl_pts = vec[0];
  ++left_ctrl_pts;

  *right_ctrl_pts = vec[last_index];
  ++right_ctrl_pts;

  while (last_index > 0)
  {
    // Construct (m - 1) control points from the m point we currently have,
    // where the new i'th point is given by: (1 - t0)*p[i] + t0*p[i + 1].
    for (i = 0; i < last_index; ++i)
    {
      vec[i] = _Point_2(comp_t0*vec[i].x() + t0*vec[i + 1].x(),
                        comp_t0*vec[i].y() + t0*vec[i + 1].y());
    }
    --last_index;

    // The first control points goes to the (front of) the right subcurve,
    // while the last control point goes to the (back of) the left subcurve.
    *left_ctrl_pts = vec[0];
    ++left_ctrl_pts;

    *right_ctrl_pts = vec[last_index];
    ++right_ctrl_pts;
  }

  // Our auxiliary vector now contains just a single entry, and this entry
  // is the curve value at t0.
  return (vec[0]);
}

} //namespace CGAL

#endif
