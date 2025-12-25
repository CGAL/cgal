// Copyright (c) 2005,2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s): Baruch Zukerman <baruchzu@post.tau.ac.il>
//            Efi Fogel       <efif@post.tau.ac.il>
//            (based on old version by Tali Zvi)

#ifndef CGAL_SURFACE_SWEEP_2_ALGORITHMS_H
#define CGAL_SURFACE_SWEEP_2_ALGORITHMS_H

#include <CGAL/license/Surface_sweep_2.h>

/*! File
 *
 * \file Definition of the surface-sweep related functions.
 */

#include <CGAL/Surface_sweep_2.h>
#include <CGAL/Do_intersect_surface_sweep_2.h>
#include <CGAL/Surface_sweep_2/Default_arr_traits.h>
#include <CGAL/Surface_sweep_2/Intersection_points_visitor.h>
#include <CGAL/Surface_sweep_2/Subcurves_visitor.h>
#include <CGAL/Surface_sweep_2/Do_interior_intersect_visitor.h>

namespace CGAL {
namespace Ss2 = Surface_sweep_2;

/*! Compute all intersection points induced by a range of input curves.
 * The intersections are calculated using the surface-sweep algorithm.
 * \param begin An input iterator for the first curve in the range.
 * \param end A input past-the-end iterator for the range.
 * \param points Output: An output iterator for the intersection points
 *                       induced by the input curves.
 * \param report_endpoints If (true), the end points of the curves are also
 *                         reported as intersection points.
 * \pre The value-type of CurveInputIterator is Traits::Curve_2, and the
 *      value-type of OutputIterator is Traits::Point_2.
 */
template <typename CurveInputIterator, typename OutputIterator, typename Traits>
OutputIterator compute_intersection_points(CurveInputIterator curves_begin,
                                           CurveInputIterator curves_end,
                                           OutputIterator points,
                                           bool report_endpoints,
                                           Traits& tr) {
  // Define the surface-sweep types:
  using Visitor = Ss2::Intersection_points_visitor<Traits, OutputIterator>;
  using Surface_sweep = Ss2::Surface_sweep_2<Visitor>;

  // Perform the sweep and obtain the intersection points.
  Visitor visitor(points, report_endpoints);
  Surface_sweep surface_sweep(&tr, &visitor);
  visitor.sweep(curves_begin, curves_end);

  return visitor.output_iterator();
}

/*!
 */
template <typename CurveInputIterator, typename OutputIterator>
OutputIterator compute_intersection_points(CurveInputIterator curves_begin,
                                           CurveInputIterator curves_end,
                                           OutputIterator points,
                                           bool report_endpoints = false) {
  using Curve = typename std::iterator_traits<CurveInputIterator>::value_type;
  typename Ss2::Default_arr_traits<Curve>::Traits   traits;
  return compute_intersection_points(curves_begin, curves_end, points, report_endpoints, traits);
}

/*! Compute all x-monotone subcurves that are disjoint in their interiors
 * induced by a range of input curves.
 * The subcurves are calculated using the surface-sweep algorithm.
 * \param begin An input iterator for the first curve in the range.
 * \param end A input past-the-end iterator for the range.
 * \param points Output: An output iterator for the subcurve.
 * \param mult_overlaps If (true), the overlapping subcurve will be reported
 *                      multiple times.
 * \pre The value-type of CurveInputIterator is Traits::Curve_2, and the
 *      value-type of OutputIterator is Traits::X_monotone_curve_2.
 */
template <typename CurveInputIterator, typename OutputIterator, typename Traits>
OutputIterator compute_subcurves(CurveInputIterator curves_begin,
                                 CurveInputIterator curves_end,
                                 OutputIterator subcurves,
                                 bool mult_overlaps, Traits& tr) {
  // Define the surface-sweep types:
  using Visitor = Ss2::Subcurves_visitor<Traits, OutputIterator>;
  using Surface_sweep = Ss2::Surface_sweep_2<Visitor>;

  // Perform the sweep and obtain the subcurves.
  Visitor visitor(subcurves, mult_overlaps);
  Surface_sweep surface_sweep(&tr, &visitor);
  visitor.sweep(curves_begin, curves_end);

  return visitor.output_iterator();
}

/*!
 */
template <typename CurveInputIterator, typename OutputIterator>
OutputIterator compute_subcurves(CurveInputIterator curves_begin,
                                 CurveInputIterator curves_end,
                                 OutputIterator subcurves,
                                 bool mult_overlaps = false) {
  using Curve = typename std::iterator_traits<CurveInputIterator>::value_type;
  typename Ss2::Default_arr_traits<Curve>::Traits traits;
  return compute_subcurves(curves_begin, curves_end, subcurves, mult_overlaps, traits);
}

/*! Determine whether any curves in a given range intersect pairwise.
 * \param begin An input iterator of the the range.
 * \param end A input past-the-end iterator of the range.
 * \return (true) if any pair of curves intersect; (false) otherwise.
 */
template <typename CurveInputIterator, typename Traits>
bool do_curves_intersect(CurveInputIterator curves_begin, CurveInputIterator curves_end, Traits& tr) {
  // Define the surface-sweep types:
  using Visitor = Ss2::Do_interior_intersect_visitor<Traits>;
  using Surface_sweep = Ss2::Do_intersect_surface_sweep_2<Visitor>;

  // Perform the sweep and obtain the subcurves.
  Visitor visitor;
  Surface_sweep surface_sweep(&tr, &visitor);
  visitor.sweep(curves_begin, curves_end);
  return visitor.found_intersection();
}

/*!
 */
template <typename CurveInputIterator>
bool do_curves_intersect(CurveInputIterator curves_begin, CurveInputIterator curves_end) {
  using Curve = typename std::iterator_traits<CurveInputIterator>::value_type;
  typename Ss2::Default_arr_traits<Curve>::Traits traits;
  return do_curves_intersect(curves_begin, curves_end, traits);
}

} // namespace CGAL

#endif
