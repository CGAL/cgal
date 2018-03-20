// Copyright (c) 2005,2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
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
#include <CGAL/Surface_sweep_2/Intersection_points_visitor.h>
#include <CGAL/Surface_sweep_2/Subcurves_visitor.h>
#include <CGAL/Surface_sweep_2/Do_interior_intersect_visitor.h>

#include <CGAL/Segment_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/Arr_circle_segment_traits_2.h>
#include <CGAL/Arr_linear_traits_2.h>

namespace CGAL {

namespace Ss2 = Surface_sweep_2;

template <typename Curve>
struct Default_arr_traits
{};

template <typename Kernel>
struct Default_arr_traits<CGAL::Segment_2<Kernel> >
{
  typedef CGAL::Arr_segment_traits_2<Kernel>                            Traits;
};

template <typename Kernel>
struct Default_arr_traits<CGAL::Arr_segment_2<Kernel> >
{
  typedef CGAL::Arr_segment_traits_2<Kernel>                            Traits;
};

template <typename SubcurveTraits>
struct Default_arr_traits<CGAL::internal::Polycurve_2
                          <SubcurveTraits, typename SubcurveTraits::Point_2> >
{
  typedef CGAL::Arr_polyline_traits_2<SubcurveTraits>                   Traits;
};

template <typename Rat_kernel_, class Alg_kernel_, class Nt_traits_>
struct Default_arr_traits<CGAL::_Conic_arc_2<Rat_kernel_, Alg_kernel_,
                                             Nt_traits_> >
{
  typedef CGAL::Arr_conic_traits_2<Rat_kernel_, Alg_kernel_, Nt_traits_>
                                                                        Traits;
};

template <typename AlgebraicKernel_d_1>
class Arr_rational_function_traits_2;

namespace Arr_rational_arc{
template <typename Algebraic_kernel_>
class Rational_arc_d_1;
}

template <typename Algebraic_kernel_>
struct Default_arr_traits<CGAL::Arr_rational_arc::
                          Rational_arc_d_1<Algebraic_kernel_> >
{
  typedef CGAL::Arr_rational_function_traits_2<Algebraic_kernel_>       Traits;
};

template <typename Kernel_, bool Filter_>
struct Default_arr_traits<CGAL::_Circle_segment_2<Kernel_, Filter_> >
{
  typedef CGAL::Arr_circle_segment_traits_2<Kernel_, Filter_>           Traits;
};

template <typename Kernel>
struct Default_arr_traits<CGAL::Arr_linear_object_2<Kernel> >
{
  typedef CGAL::Arr_linear_traits_2<Kernel>                             Traits;
};

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
                                           Traits &tr)
{
  // Define the surface-sweep types:
  typedef Ss2::Intersection_points_visitor<Traits, OutputIterator>
                                                                Visitor;
  typedef Ss2::Surface_sweep_2<Visitor>                         Surface_sweep;

  // Perform the sweep and obtain the intersection points.
  Visitor visitor(points, report_endpoints);
  Surface_sweep surface_sweep(&tr, &visitor);
  visitor.sweep(curves_begin, curves_end);

  return visitor.output_iterator();
}

template <typename CurveInputIterator, typename OutputIterator>
OutputIterator compute_intersection_points(CurveInputIterator curves_begin,
                                           CurveInputIterator curves_end,
                                           OutputIterator points,
                                           bool report_endpoints = false)
{
  typedef typename std::iterator_traits<CurveInputIterator>::value_type  Curve;

  typename Default_arr_traits<Curve>::Traits   traits;

  return compute_intersection_points(curves_begin, curves_end, points,
                                     report_endpoints, traits);
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
                                 bool mult_overlaps, Traits& tr)
{
  // Define the surface-sweep types:
  typedef Ss2::Subcurves_visitor<Traits, OutputIterator>        Visitor;
  typedef Ss2::Surface_sweep_2<Visitor>                         Surface_sweep;

  // Perform the sweep and obtain the subcurves.
  Visitor visitor(subcurves, mult_overlaps);
  Surface_sweep surface_sweep(&tr, &visitor);
  visitor.sweep(curves_begin, curves_end);

  return visitor.output_iterator();
}

template <typename CurveInputIterator, typename OutputIterator>
OutputIterator compute_subcurves(CurveInputIterator curves_begin,
                                 CurveInputIterator curves_end,
                                 OutputIterator subcurves,
                                 bool mult_overlaps = false)
{
  typedef typename std::iterator_traits<CurveInputIterator>::value_type  Curve;
  typename Default_arr_traits<Curve>::Traits m_traits;
  return compute_subcurves(curves_begin, curves_end, subcurves, mult_overlaps,
                           m_traits);
}

/*! Determine if there occurs an intersection between any pair of curves in
 * a given range.
 * \param begin An input iterator for the first curve in the range.
 * \param end A input past-the-end iterator for the range.
 * \return (true) if any pair of curves intersect; (false) otherwise.
 */
template <typename CurveInputIterator, typename Traits>
bool do_curves_intersect(CurveInputIterator curves_begin,
                         CurveInputIterator curves_end, Traits& tr)
{
  // Define the surface-sweep types:
  typedef Ss2::Do_interior_intersect_visitor<Traits>            Visitor;
  typedef Ss2::Surface_sweep_2<Visitor>                         Surface_sweep;

  // Perform the sweep and obtain the subcurves.
  Visitor visitor;
  Surface_sweep surface_sweep(&tr, &visitor);
  visitor.sweep(curves_begin, curves_end);
  return visitor.found_intersection();
}

template <typename CurveInputIterator>
bool do_curves_intersect(CurveInputIterator curves_begin,
                         CurveInputIterator curves_end)
{
  typedef typename std::iterator_traits<CurveInputIterator>::value_type  Curve;

  typename Default_arr_traits<Curve>::Traits m_traits;
  return do_curves_intersect(curves_begin, curves_end, m_traits);
}

} // namespace CGAL

#endif
