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
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Ron Wein        <wein@post.tau.ac.il>
//                 Efi Fogel       <efif@post.tau.ac.il>

#ifndef CGAL_SURFACE_SWEEP_2_UTILS_H
#define CGAL_SURFACE_SWEEP_2_UTILS_H

#include <CGAL/license/Surface_sweep_2.h>

/*! \file
 *
 * Auxiliary functions for the usage of the various sweep-line visitors.
 */

#include <CGAL/basic.h>
#include <CGAL/Object.h>
#include <CGAL/assertions.h>
#include <vector>
#include <algorithm>
#include <CGAL/Arr_enums.h>

namespace CGAL {
namespace Surface_sweep_2 {

/*! Subdivide a range of input curves into x-monotone objects.
 * \param begin The first input curve (of type Curve_2).
 * \param end A part-the-end iterator for the input curves.
 * \param x_curves Output: The x-monotone subcurves
 *                         (of type X_monotone_curve_2).
 * \param iso_points Output: The isolated points (of type Point_2).
 * \param tr A geometry-traits class.
 */
template <typename Traits,
          typename CurveInputIter,
          typename XCurveOutIter,
          typename PointOutIter>
void make_x_monotone(CurveInputIter begin, CurveInputIter end,
                     XCurveOutIter x_curves,
                     PointOutIter iso_points,
                     const Traits* tr)
{
  // Split the input curves into x-monotone objects.
  std::size_t num_of_curves = std::distance(begin, end);
  std::vector<Object> object_vec;
  CurveInputIter iter;

  object_vec.reserve(num_of_curves);
  for (iter = begin; iter != end; ++iter) {
    tr->make_x_monotone_2_object()(*iter, std::back_inserter(object_vec));
  }

  // Transform each object to either a point or an x-monotone curve.
  typedef typename Traits::X_monotone_curve_2    X_monotone_curve_2;
  typedef typename Traits::Point_2               Point_2;

  const X_monotone_curve_2* xcv;
  const Point_2* pt;
  unsigned int i;

  for (i = 0 ; i < object_vec.size() ; ++i) {
    xcv = object_cast<X_monotone_curve_2> (&(object_vec[i]));

    if (xcv != NULL) {
      // The object is an x-monotone curve.
      *x_curves = *xcv;
      ++x_curves;
    }
    else {
      // The object is an isolated point.
      pt = object_cast<Point_2> (&(object_vec[i]));
      CGAL_assertion (pt != NULL);

      *iso_points = *pt;
      ++iso_points;
    }
  }
}

/*! Given an arrangement and two ranges of x-monotone curves and isolated
 * points, representing objects that should be inserted into the arrangement,
 * create two output sets of extended x-monotone curves and isolated points,
 * including the arrangement edges and isolated vertices.
 * \param arr The input arrangement.
 * \param xcvs_begin The first input x-monotone curve
 *                   (of type Arrangement::X_monotone_Curve_2).
 * \param xcvs_end A past-the-end iterator for the input x-monotone curves.
 * \param pts_begin The first isolated input point
 *                   (of type Arrangement::Point_2).
 * \param pts_end A past-the-end iterator for the isolated input points.
 * \param x_curves Output: The extended x-monotone subcurves
 *                         (of type ExTraits::X_monotone_curve_2).
 * \param iso_points Output: The extended isolated points
 *                           (of type ExTraits::Point_2).
 * \param ex_tr An extended geometry-traits class.
 *              This parameter is not actually in use, but is needed in order
 *              to instantiate the template parameter ExTraits.
 */
template <typename Arrangement,
          typename ExTraits,
          typename XCurveInputIter,
          typename PointInputIter,
          typename XCurveOutIter,
          typename PointOutIter>
void prepare_for_sweep(Arrangement& arr,
                       XCurveInputIter xcvs_begin, XCurveInputIter xcvs_end,
                       PointInputIter pts_begin, PointInputIter pts_end,
                       XCurveOutIter x_curves,
                       PointOutIter iso_points,
                       const ExTraits * /* ex_tr */)
{
  typedef typename Arrangement::Vertex_iterator       Vertex_iterator;
  typedef typename Arrangement::Edge_iterator         Edge_iterator;
  typedef typename Arrangement::Vertex_handle         Vertex_handle;
  typedef typename Arrangement::Halfedge_handle       Halfedge_handle;

  typedef typename ExTraits::X_monotone_curve_2       Ex_x_monotone_curve_2;
  typedef typename ExTraits::Point_2                  Ex_point_2;

  // Go over the input objects and copy them to the output iterators.
  XCurveInputIter xcv_it;
  PointInputIter pt_it;

  for (xcv_it = xcvs_begin; xcv_it != xcvs_end; ++xcv_it) {
    *x_curves = Ex_x_monotone_curve_2 (*xcv_it);
    ++x_curves;
  }

  for (pt_it = pts_begin; pt_it != pts_end; ++pt_it) {
    *iso_points = Ex_point_2 (*pt_it);
    ++iso_points;
  }

  // Go over the arrangement edges and insert their associated x-monotone
  // curves into the output iterator. To each curve we attach a handle to the
  // halfedge that goes from right to left.
  Edge_iterator     eit;
  Halfedge_handle   he;

  for (eit = arr.edges_begin(); eit != arr.edges_end(); ++eit) {
    if (eit->direction() == ARR_LEFT_TO_RIGHT) he = eit->twin();
    else he = eit;

    *x_curves = Ex_x_monotone_curve_2 (he->curve(), he);
    ++x_curves;
  }

  // Go over the isolated arrangement vertices and insert their associated
  // points into the output iterator. To each point we attach a handle to its
  // vertex.
  Vertex_iterator   vit;
  Vertex_handle     v;

   for (vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit) {
     v = vit;
     if (v->is_isolated()) {
       *iso_points = Ex_point_2 (v->point(), v);
       ++iso_points;
     }
   }
}

} // namespace CGAL
} // namespace Surface_sweep_2

#endif
