// Copyright (c) 2006,2007,2009,2010,2011,2026 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s) : Baruch Zukerman <baruchzu@post.tau.ac.il>
//             Ron Wein        <wein@post.tau.ac.il>
//             Efi Fogel       <efif@post.tau.ac.il>

#ifndef CGAL_SURFACE_SWEEP_2_UTILS_H
#define CGAL_SURFACE_SWEEP_2_UTILS_H

#include <CGAL/license/Surface_sweep_2.h>

/*! \file
 *
 * Auxiliary functions for the usage of the various sweep-line visitors.
 */

#include <vector>
#include <algorithm>

#include <boost/iterator/function_output_iterator.hpp>

#include <CGAL/basic.h>
#include <CGAL/assertions.h>
#include <CGAL/Arr_enums.h>
#include <CGAL/variant_output_iterator.h>

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
template <typename CurveInputIterator, typename XCurveOutputIterator, typename PointOutputIterator, typename Traits>
void make_x_monotone(CurveInputIterator begin, CurveInputIterator end,
                     XCurveOutputIterator it_xcv, PointOutputIterator it_pt, const Traits& traits) {
  using Point_2 = typename Traits::Point_2;
  using X_monotone_curve_2 = typename Traits::X_monotone_curve_2;
  using Variant = std::variant<Point_2, X_monotone_curve_2>;
  auto out =
    CGAL::make_variant_output_iterator<Variant>([it_pt](const Point_2& pt) mutable { *it_pt++ = pt; },
                                                [it_xcv](const X_monotone_curve_2& xcv) mutable { *it_xcv++ = xcv; });
  auto mk_x_monotone = traits.make_x_monotone_2_object();
  for (auto it = begin; it != end; ++it) mk_x_monotone(*it, out);
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
 * \param xcurves Output: The extended x-monotone subcurves
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
                       XCurveOutIter xcurves, PointOutIter iso_points, const ExTraits * /* ex_traits */) {
  using Vertex_handle = typename Arrangement::Vertex_handle;
  using Halfedge_handle = typename Arrangement::Halfedge_handle;

  using Ex_x_monotone_curve_2 = typename ExTraits::X_monotone_curve_2;
  using Ex_point_2 = typename ExTraits::Point_2;

  // Go over the input objects and copy them to the output iterators.
  for (auto xcv_it = xcvs_begin; xcv_it != xcvs_end; ++xcv_it) *xcurves++ = Ex_x_monotone_curve_2(*xcv_it);
  for (auto pt_it = pts_begin; pt_it != pts_end; ++pt_it) *iso_points++ = Ex_point_2(*pt_it);

  // Go over the arrangement edges and insert their associated x-monotone
  // curves into the output iterator. To each curve we attach a handle to the
  // halfedge that goes from right to left.
  for (auto eit = arr.edges_begin(); eit != arr.edges_end(); ++eit) {
    Halfedge_handle he = (eit->direction() == ARR_LEFT_TO_RIGHT) ? eit->twin() : eit;
    *xcurves++ = Ex_x_monotone_curve_2(he->curve(), he);
  }

  // Go over the isolated arrangement vertices and insert their associated
  // points into the output iterator. To each point we attach a handle to its
  // vertex.
  for (auto vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit) {
    Vertex_handle v = vit;
    if (v->is_isolated()) *iso_points++ = Ex_point_2(v->point(), v);
  }
}

} // namespace CGAL
} // namespace Surface_sweep_2

#endif
