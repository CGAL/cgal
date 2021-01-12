// Copyright (c) 2007,2008,2009,2010,2011 Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Eric Berberich <eric@mpi-inf.mpg.de>
//                 Arno Eigenwillig <arno@mpi-inf.mpg.de>
//                 Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

#ifndef CGAL_CURVED_KERNEL_VIA_ANALYSIS_2_MAKE_X_MONOTONE_2_H
#define CGAL_CURVED_KERNEL_VIA_ANALYSIS_2_MAKE_X_MONOTONE_2_H

/*!\file include/CGAL/Curved_kernel_via_analysis_2/Make_x_monotone_2.h
 * \brief Defines \c Make_x_monotone_2 functor
 */

#include <CGAL/config.h>
#include <CGAL/iterator.h>
#include <CGAL/Handle_with_policy.h>

// TODO remove polynomial_traits
#include <CGAL/Polynomial_traits_d.h>

namespace CGAL {

namespace internal {

/*!\brief
 * Splits a curve that can be analyzed
 * into connected x-monotone sweepable arcs and isolated points.
 *
 * Arcs are stored as CurvedKernelViaAnalysis_2::Arc_2 objects, and
 * each is either vertical or consists of an x-monotone piece
 * of constant arc number wrt to the curve at every interior x-coordinate.
 * Isolated points are stored as \c CurvedKernelViaAnalysis_2::Point_2 objects.
 *
 * The resulting arcs and points are written to the output iterator as
 * polymorphic \c variant. Past-the-end value of the iterator is returned.
 *
 * EF: I believe that the inheritance from binary_function is not exploited,
 *     and thus redundant, but I keep it anyway.
 */
template < class CurvedKernelViaAnalysis_2,
           class ConstructArc_2 =
           typename CurvedKernelViaAnalysis_2::Construct_arc_2 >
struct Make_x_monotone_2 :
    public CGAL::cpp98::binary_function<
        typename CurvedKernelViaAnalysis_2::Curve_2,
        CGAL::cpp98::iterator<std::output_iterator_tag,
                              boost::variant<
                                typename CurvedKernelViaAnalysis_2::Point_2,
                                typename CurvedKernelViaAnalysis_2::Arc_2> >,
        CGAL::cpp98::iterator<std::output_iterator_tag,
                              boost::variant<
                                typename CurvedKernelViaAnalysis_2::Point_2,
                                typename CurvedKernelViaAnalysis_2::Arc_2> > >
{

    //!\name Public types
    //!@{

    //! this instance's first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! this instance's second template parameter
    typedef ConstructArc_2 Construct_arc_2;

    //! type of curve kernel
    typedef typename Curved_kernel_via_analysis_2::Curve_kernel_2
    Curve_kernel_2;

    //! type of x-coordinate
    typedef typename Curve_kernel_2::Coordinate_1 Coordinate_1;

    //! type of xy-coordinate
    typedef typename Curve_kernel_2::Coordinate_2 Coordinate_2;

    //! type of curve analysis
    typedef typename Curve_kernel_2::Curve_analysis_2 Curve_analysis_2;

    //! type of vertical line
    typedef typename Curve_analysis_2::Status_line_1 Status_line_1;

    //! type of point on curve
    typedef typename Curved_kernel_via_analysis_2::Point_2 Point_2;

    //! type of curve arc
    typedef typename Curved_kernel_via_analysis_2::Arc_2 Arc_2;

    //! type of not necessarily x-monotone curve arc
    typedef typename Curved_kernel_via_analysis_2::Non_x_monotone_arc_2
        Non_x_monotone_arc_2;

  //!@}

  //!\name Constructors
  //!@{

  /*!\brief
   * Standard constructor
   *
   * \param kernel The kernel instance to use
   */
  Make_x_monotone_2(Curved_kernel_via_analysis_2 *kernel) :
    _m_curved_kernel(kernel)
  {
    CGAL_assertion(kernel != nullptr);
  }

  //!@}

  //!\name Functor invokation
  //!@{

  // TODO add operator for non-x-monotone arc

  /*!\brief
   * Splits a curve into x-monotone arcs and isolated points
   *
   * \param curve The input curve
   * \param oi the output iterator for the result. Its dereference type is a
   *           variant that wraps a Point_2 or an X_monotone_curve_2 objects.
   * \return Past-the-end iterator of \c oi
   */
  template <typename OutputIterator>
  OutputIterator operator()(Curve_analysis_2 curve, OutputIterator oi)
  {
    typedef boost::variant<Point_2, Arc_2>      Make_x_monotone_result;

    Construct_arc_2 construct_arc_2 =
      _m_curved_kernel->construct_arc_2_object();
    // use CGAL::Total_degree ?
    if (typename CGAL::Polynomial_traits_d<
        typename Curve_analysis_2::Polynomial_2 >::
        Total_degree()(curve.polynomial_2()) < 1)
    {
      return oi;
    }

    Status_line_1 evt_line1, evt_line2,
      int_line = curve.status_line_of_interval(0);
    int total_events = curve.number_of_status_lines_with_event();
    // handle special case of a curve without any events
    if(total_events == 0) {
      for (int k = 0; k < int_line.number_of_events(); k++)
        *oi++ = Make_x_monotone_result(construct_arc_2(curve, k));
      return oi;
    }
    _m_curve = curve;
    typedef typename Curved_kernel_via_analysis_2::
      Curve_interval_arcno_cache CIA_cache;
    const CIA_cache& map_interval_arcno =
      _m_curved_kernel->interval_arcno_cache();

    typename Curved_kernel_via_analysis_2::Construct_point_2
      construct_point =
      _m_curved_kernel->construct_point_2_object();

    typename CIA_cache::result_type info1, info2;
    std::vector<Point_2> min_pts, max_pts;
    Coordinate_1 min_x, max_x;
    int i, k, n;
    Arc_2 arc;
    // first handle segments before first event
    evt_line1 = curve.status_line_at_event(0);
    max_x = evt_line1.x();

    for (k = 0; k < evt_line1.number_of_events(); k++)
      max_pts.push_back(construct_point(max_x, curve, k));

    //std::cout << "handling events over the 1st interval\n";
    for (k = 0; k < int_line.number_of_events(); k++) {

      info1 = map_interval_arcno(evt_line1, 1, k);
      if (info1.second != CGAL::ARR_INTERIOR) {
        arc = construct_arc_2(CGAL::ARR_MIN_END, max_x,
                              (info1.second ==
                               CGAL::ARR_BOTTOM_BOUNDARY ?
                               CGAL::ARR_MIN_END : CGAL::ARR_MAX_END),
                              curve, k);
      }
      else {
        arc = construct_arc_2(max_pts[info1.first], CGAL::ARR_MIN_END,
                              curve, k, info1.first);
      }
      *oi++ = Make_x_monotone_result(arc);
    }
    min_pts = max_pts;
    max_pts.clear();
    min_x = max_x;

    // next handle arcs between events, including isolated points
    for (i = 0; i < total_events-1; i++) {
      evt_line1 = curve.status_line_at_event(i);
      evt_line2 = curve.status_line_at_event(i+1);
      max_x = evt_line2.x();
      oi = _handle_vertical_and_isolated(evt_line1, min_x, min_pts, oi);

      n = evt_line2.number_of_events();
      for (k = 0; k < n; k++)
        max_pts.push_back(construct_point(max_x, curve, k));

      n = curve.status_line_of_interval(i+1).number_of_events();
      CGAL::Arr_curve_end inf1_end, inf2_end;
      for (k = 0; k < n; k++) {
        info1 = map_interval_arcno(evt_line1, 0, k);
        info2 = map_interval_arcno(evt_line2, 1, k);
        inf2_end = (info2.second == CGAL::ARR_BOTTOM_BOUNDARY ?
                    CGAL::ARR_MIN_END : CGAL::ARR_MAX_END);

        if (info1.second != CGAL::ARR_INTERIOR) {
          inf1_end = (info1.second == CGAL::ARR_BOTTOM_BOUNDARY ?
                      CGAL::ARR_MIN_END : CGAL::ARR_MAX_END);
          if (info2.second != CGAL::ARR_INTERIOR) {
            arc = construct_arc_2(min_x, inf1_end, max_x, inf2_end,
                                  curve, k);
          }
          else {
            arc = construct_arc_2(max_pts[info2.first], min_x,
                                  inf1_end, curve, k, info2.first);
          }
        }
        else if (info2.second != CGAL::ARR_INTERIOR) {
          arc = construct_arc_2(min_pts[info1.first],  max_x,
                                inf2_end, curve, k, info1.first);
        }
        else {
          arc = construct_arc_2(min_pts[info1.first],
                                max_pts[info2.first],
                                curve, k, info1.first, info2.first);
        }
        *oi++ = Make_x_monotone_result(arc);
      }
      min_pts = max_pts;
      max_pts.clear();
      min_x = max_x;
    }

    // here: min_x/min_pts hold information about the last event line
    // event_line2 - points to the last event line
    // vertical line or isolated points at last event?
    evt_line2 = curve.status_line_at_event(total_events-1);
    min_x = evt_line2.x();
    oi = _handle_vertical_and_isolated(evt_line2, min_x, min_pts, oi);

    n = curve.status_line_of_interval(total_events).number_of_events();
    for (k = 0; k < n; k++) {
      info1 = map_interval_arcno(evt_line2, 0, k);
      if (info1.second != CGAL::ARR_INTERIOR) {
        arc = construct_arc_2(CGAL::ARR_MAX_END, min_x,
                              (info1.second == CGAL::ARR_BOTTOM_BOUNDARY ?
                               CGAL::ARR_MIN_END : CGAL::ARR_MAX_END), curve, k
                              );
      }
      else {
        arc = construct_arc_2(min_pts[info1.first],
                              CGAL::ARR_MAX_END, curve, k,
                              info1.first);
      }
      *oi++ = Make_x_monotone_result(arc);
    }
    return oi;
  }
  //!@}

private:
  //!\name Private members
  //!@{

  /*!\brief
   * Constructs vertical arcs and isolated points at event line
   *
   * \param cv_line The event line in focus
   * \param x x-coordinate of event
   * \param pts Points at event line
   * \param oi the output iterator for the result. Its dereference type is a
   *           variant that wraps a Point_2 or an X_monotone_curve_2 objects.
   * \return Past-the-end iterator of \c oi
   */
  template <class OutputIterator>
  OutputIterator
  _handle_vertical_and_isolated(Status_line_1 cv_line,
                                Coordinate_1 x, std::vector<Point_2> pts,
                                OutputIterator oi) const
  {
    typedef boost::variant<Point_2, Arc_2>      Make_x_monotone_result;

    Construct_arc_2 construct_arc_2 =
      _m_curved_kernel->construct_arc_2_object();

    int n = cv_line.number_of_events(), j;
    if (cv_line.covers_line()) { // look for vertical arcs
      if (n > 0) {
        // the first vertical ray
        *oi++ =
          Make_x_monotone_result(construct_arc_2(pts[0], CGAL::ARR_MIN_END,
                                                 _m_curve));
        for (j = 0; j < n-1; j++)  // interior bounded arcs
          *oi++ = Make_x_monotone_result(construct_arc_2(pts[j], pts[j+1],
                                                    _m_curve));
        // the last vertical ray
        *oi++ =
          Make_x_monotone_result(construct_arc_2(pts[n-1], CGAL::ARR_MAX_END,
                                                 _m_curve));
      }
      else // unbounded vertical line
        *oi++ = Make_x_monotone_result(construct_arc_2(x, _m_curve));
      return oi;
    }
    // look for isolated points
    std::pair<int, int> ipair;
    for (j = 0; j < n; j++) {
      ipair = cv_line.number_of_incident_branches(j);
      if (ipair.first == 0&&ipair.second == 0) {
        //std::cout << "isolated point found\n";
        typename Curved_kernel_via_analysis_2::Construct_point_2
          construct_point = _m_curved_kernel->construct_point_2_object();

        *oi++ = Make_x_monotone_result(construct_point(x, _m_curve, j));
      }
    }
    return oi;
  }

  //!@}

  //!\name Private data
  //!@{

  //! pointer to \c Curved_kernel_via_analysis_2
  Curved_kernel_via_analysis_2 *_m_curved_kernel;

  //! to avoid passing curve as a parameter
  Curve_analysis_2 _m_curve;

  //!@}
}; // struct Make_x_monotone

} // namespace internal

} //namespace CGAL

#endif // CGAL_CURVED_KERNEL_VIA_ANALYSIS_2_MAKE_X_MONOTONE_2_H
//EOF
