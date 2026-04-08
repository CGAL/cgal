// Copyright (c) 2005,2006,2007,2009,2010,2011,2026 Tel-Aviv University (Israel).
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
//             Efi Fogel       <efif@post.tau.ac.il>
//             (based on old version by Tali Zvi)

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
#include <CGAL/Surface_sweep_2/Do_intersect_visitor.h>
#include <CGAL/Surface_sweep_2/Do_interior_intersect_visitor.h>
#include <CGAL/Surface_sweep_2/Intersection_polylines_visitor.h>
#include <CGAL/Surface_sweep_2/Surface_sweep_2_utils.h>

#include <CGAL/Arr_curve_data_traits_2.h>

namespace CGAL {
namespace Surface_sweep_2 {

/*! Determine whether any curves in a given range intersect pairwise.
 * \param begin An input iterator of the the range.
 * \param end A input past-the-end iterator of the range.
 * \param consider_common_endpoints Indicates whether common endpoints should be considered.
 * \return (true) if any pair of curves intersect; (false) otherwise.
 */
template <typename CurveInputIterator, typename Traits>
bool do_intersect(CurveInputIterator begin, CurveInputIterator end, bool consider_common_endpoints, Traits& traits) {
  // If the curves are not \f$x\f$-monotone, subdivide them into \f$x\f$-monotone curves.
  using Visitor = Do_intersect_visitor<Traits>;
  using Surface_sweep = Do_intersect_surface_sweep_2<Visitor>;

  Visitor visitor;
  Surface_sweep surface_sweep(&traits, &visitor, consider_common_endpoints);

  using X_monotone_curve_2 = typename Traits::X_monotone_curve_2;
  using value_type = typename std::iterator_traits<CurveInputIterator>::value_type;
  if constexpr (std::is_same_v<value_type, X_monotone_curve_2>) {
    surface_sweep.do_intersect_sweep(begin, end);
    return visitor.do_intersect();
  }
  else {
    using Point_2 = typename Traits::Point_2;
    std::vector<X_monotone_curve_2> xcurves;
    std::vector<Point_2> points;
    xcurves.reserve(std::distance(begin, end));
    make_x_monotone(begin, end, std::back_inserter(xcurves), std::back_inserter(points), traits);
    surface_sweep.do_intersect_sweep(xcurves.begin(), xcurves.end(), points.begin(), points.end());
    return visitor.do_intersect();
  }
}

/*!
 */
template <typename CurveInputIterator>
bool do_intersect(CurveInputIterator begin, CurveInputIterator end, bool consider_common_endpoints = true) {
  using Curve = typename std::iterator_traits<CurveInputIterator>::value_type;
  typename Default_arr_traits<Curve>::Traits traits;
  return do_intersect(begin, end, consider_common_endpoints, traits);
}

} // namespace Surface_sweep_2

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
OutputIterator compute_intersection_points(CurveInputIterator begin, CurveInputIterator end,
                                           OutputIterator points, bool report_endpoints, Traits& traits) {
  // Define the surface-sweep types:
  using Visitor = Ss2::Intersection_points_visitor<Traits, OutputIterator>;
  using Surface_sweep = Ss2::Surface_sweep_2<Visitor>;

  // Perform the sweep and obtain the intersection points.
  Visitor visitor(points, report_endpoints);
  Surface_sweep surface_sweep(&traits, &visitor);

  using X_monotone_curve_2 = typename Traits::X_monotone_curve_2;
  using value_type = typename std::iterator_traits<CurveInputIterator>::value_type;
  if constexpr (std::is_same_v<value_type, X_monotone_curve_2>) {
    surface_sweep.sweep(begin, end);
    return visitor.output_iterator();
  }
  else {
    using Point_2 = typename Traits::Point_2;
    std::vector<X_monotone_curve_2> xcurves;
    std::vector<Point_2> points;
    xcurves.reserve(std::distance(begin, end));
    Ss2::make_x_monotone(begin, end, std::back_inserter(xcurves), std::back_inserter(points), traits);
    surface_sweep.sweep(xcurves.begin(), xcurves.end(), points.begin(), points.end());
    return visitor.output_iterator();
  }
}

/*!
 */
template <typename CurveInputIterator, typename OutputIterator>
OutputIterator compute_intersection_points(CurveInputIterator begin, CurveInputIterator end,
                                           OutputIterator points, bool report_endpoints = false) {
  using Curve = typename std::iterator_traits<CurveInputIterator>::value_type;
  typename Ss2::Default_arr_traits<Curve>::Traits traits;
  return compute_intersection_points(begin, end, points, report_endpoints, traits);
}

/*! Compute all x-monotone subcurves that are disjoint in their interiors
 * induced by a range of input curves.
 * The subcurves are calculated using the surface-sweep algorithm.
 * \param begin An input iterator for the first curve in the range.
 * \param end A input past-the-end iterator for the range.
 * \param subcurves Output: An output iterator for the subcurve.
 * \param mult_overlaps If (true), the overlapping subcurve will be reported
 *                      multiple times.
 * \pre The value-type of `CurveInputIterator` is `Traits::Curve_2`, and the
 *      value-type of `OutputIterator` is `Traits::X_monotone_curve_2`.
 */
template <typename CurveInputIterator, typename OutputIterator, typename Traits>
OutputIterator compute_subcurves(CurveInputIterator begin, CurveInputIterator end,
                                 OutputIterator subcurves, bool mult_overlaps, Traits& traits) {
  // Define the surface-sweep types:
  using Visitor = Ss2::Subcurves_visitor<Traits, OutputIterator>;
  using Surface_sweep = Ss2::Surface_sweep_2<Visitor>;

  // Perform the sweep and obtain the subcurves.
  Visitor visitor(subcurves, mult_overlaps);
  Surface_sweep surface_sweep(&traits, &visitor);

  using X_monotone_curve_2 = typename Traits::X_monotone_curve_2;
  using value_type = typename std::iterator_traits<CurveInputIterator>::value_type;
  if constexpr (std::is_same_v<value_type, X_monotone_curve_2>) {
    surface_sweep.sweep(begin, end);
    return visitor.output_iterator();
  }
  else {
    using Point_2 = typename Traits::Point_2;
    std::vector<X_monotone_curve_2> xcurves;
    std::vector<Point_2> points;
    xcurves.reserve(std::distance(begin, end));
    Ss2::make_x_monotone(begin, end, std::back_inserter(xcurves), std::back_inserter(points), traits);
    surface_sweep.sweep(xcurves.begin(), xcurves.end(), points.begin(), points.end());
    return visitor.output_iterator();
  }
}

/*!
 */
template <typename CurveInputIterator, typename OutputIterator>
OutputIterator compute_subcurves(CurveInputIterator begin, CurveInputIterator end,
                                 OutputIterator subcurves, bool mult_overlaps = false) {
  using Curve = typename std::iterator_traits<CurveInputIterator>::value_type;
  typename Ss2::Default_arr_traits<Curve>::Traits traits;
  return compute_subcurves(begin, end, subcurves, mult_overlaps, traits);
}

/*! Determine whether any curves in a given range intersect pairwise.
 * \param begin An input iterator of the the range.
 * \param end A input past-the-end iterator of the range.
 * \return (true) if any pair of curves intersect; (false) otherwise.
 */
template <typename CurveInputIterator, typename Traits>
CGAL_DEPRECATED
bool do_curves_intersect(CurveInputIterator begin, CurveInputIterator end, Traits& traits)
{ return Ss2::do_intersect(begin, end, false, traits); }

/*!
 */
template <typename CurveInputIterator>
CGAL_DEPRECATED
bool do_curves_intersect(CurveInputIterator begin, CurveInputIterator end)
{ return Ss2::do_intersect(begin, end, false); }

/*! Subdivide a range of input curves according to their pairwise intersections.
 * Each curve is subdivided into sub-curves, referred to as polylines.
 *
 * \tparam CurveInputIterator model of `ForwardIterator` whose `value_type` is `Traits::Curve_2`
 * \tparam PointRange model of `RandomAccessContainer` and `BackInsertionSequence` whose value type is `Traits::Point_2`
 * \tparam PolylineRange model of `RandomAccessContainer` and `BackInsertionSequence` whose `value_type` is itself a model of `RandomAccessConatiner` and `BackInsertionSequence`
 *                       whose `value_type` is an unsigned integer type convertible to `std::size_t`.
 *
 * \param begin An input iterator for the first curve in the range.
 * \param end A input past-the-end iterator for the range.
 * \param output_points A range that will be populated with all vertices induced by the input curves (their endpoints and intersection points).
 * \param output_polylines A range that will be populated with index sequence referring to
 *   `output_points`. Each sequence corresponds to one input curve and
 *   describes the polyline obtained after subdividing the curve at
 *   intersection points.
 * \pre The value-type of `CurveInputIterator` is `Traits::Curve_2`, and the
 *      value-type of `OutputIterator` is `Traits::X_monotone_curve_2`.
 */
template <typename CurveInputIterator, typename PointRange, typename PolylineRange, typename Traits>
void compute_intersection_polylines(CurveInputIterator begin,
                                    CurveInputIterator end,
                                    PointRange& output_points,
                                    PolylineRange& output_polylines,
                                    Traits& traits)
{
  using Internal_traits = Arr_curve_data_traits_2<Traits, std::size_t>;
  using Visitor = Ss2::Intersection_polylines_visitor<Internal_traits, PointRange, PolylineRange>;
  using Surface_sweep = Ss2::Surface_sweep_2<Visitor>;

  using Internal_curve_2 = typename Internal_traits::Curve_2;

  using X_monotone_curve_2 = typename Traits::X_monotone_curve_2;
  using value_type = typename std::iterator_traits<CurveInputIterator>::value_type;

  std::vector<Internal_curve_2> curves;
  curves.reserve(std::distance(begin, end));
  std::size_t i=0;
  for(auto it=begin; it!=end; ++it)
    curves.emplace_back(*it, i++);

  Internal_traits internal_traits(traits);
  output_points.reserve(2 * std::distance(begin, end));
  output_polylines.resize(std::distance(begin, end));
  Visitor visitor(output_points, output_polylines);
  Surface_sweep surface_sweep(&internal_traits, &visitor);

  if constexpr (std::is_same_v<value_type, X_monotone_curve_2>) {
    surface_sweep.sweep(curves.begin(), curves.end());
  } else {
    using Internal_X_curve_2 = typename Internal_traits::X_monotone_curve_2;
    using Point_2 = typename Internal_traits::Point_2;

    std::vector<Internal_X_curve_2> xcurves;
    std::vector<Point_2> points;
    xcurves.reserve(std::distance(begin, end));
    Ss2::make_x_monotone(curves.begin(), curves.end(), std::back_inserter(xcurves), std::back_inserter(points), internal_traits);
    surface_sweep.sweep(xcurves.begin(), xcurves.end(), points.begin(), points.end());
  }
}

template <typename CurveInputIterator, typename PointRange, typename PolylineRange>
void compute_intersection_polylines(CurveInputIterator begin,
                                    CurveInputIterator end,
                                    PointRange& output_points,
                                    PolylineRange& output_polylines)
{
  typedef typename std::iterator_traits<CurveInputIterator>::value_type  Curve;

  typename Ss2::Default_arr_traits<Curve>::Traits m_traits;
  compute_intersection_polylines(begin, end, output_points, output_polylines, m_traits);
}

} // namespace CGAL

#endif
