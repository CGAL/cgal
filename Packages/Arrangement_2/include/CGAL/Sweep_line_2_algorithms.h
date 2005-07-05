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
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 (based on old version by Tali Zvi)
#ifndef CGAL_SWEEP_LINE_2_ALGORITHMS_H
#define CGAL_SWEEP_LINE_2_ALGORITHMS_H

/*!
 * \file Definition of the sweep-line related functions.
 */

#include <CGAL/Sweep_line_2.h>
#include <CGAL/Sweep_line_2/Sweep_line_points_visitor.h>
#include <CGAL/Sweep_line_2/Sweep_line_subcurves_visitor.h>
#include <CGAL/Sweep_line_2/Sweep_line_do_curves_x_visitor.h>

CGAL_BEGIN_NAMESPACE

/*!
 * Compute all intersection points induced by a range of input curves.
 * The intersections are calculated using the sweep-line algorithm.
 * \param begin An input iterator for the first curve in the range.
 * \param end A input past-the-end iterator for the range.
 * \param points Output: An output iterator for the intersection points 
 *                       induced by the input curves.
 * \param report_endpoints If (true), the end points of the curves are also
 *                         reported as intersection points.
 *                         Defaults value is (false).
 * \pre The value-type of CurveInputIterator is Traits::Curve_2, and the
 *      value-type of OutputIterator is Traits::Point_2.
 */
template <class Traits, class CurveInputIterator, class OutputIterator>
OutputIterator get_intersection_points (CurveInputIterator curves_begin,
                                        CurveInputIterator curves_end,
                                        OutputIterator points,
                                        Traits &tr,
                                        bool report_endpoints = false)
{
  // Define the sweep-line types:
  typedef Sweep_line_points_visitor<Traits,OutputIterator>  Visitor;
  typedef Sweep_line_2< Traits,
                        Sweep_line_points_visitor<Traits,
                                                  OutputIterator> >
                                                            Sweep_line;

  // Perform the sweep and obtain the intersection points.
  Visitor     visitor (points, report_endpoints, &tr);
  Sweep_line  sweep_line (&tr, &visitor);
  visitor.sweep(curves_begin, curves_end);

  return (visitor.get_output_iterator());
}

/*!
 * Compute all x-monotone subcurves that are disjoint in their interiors
 * induced by a range of input curves.
 * The subcurves are calculated using the sweep-line algorithm.
 * \param begin An input iterator for the first curve in the range.
 * \param end A input past-the-end iterator for the range.
 * \param points Output: An output iterator for the subcurve. 
 * \param mult_overlaps If (true), the overlapping subcurve will be reported
 *                      multiple times.
 *                      Defaults value is (false).
 * \pre The value-type of CurveInputIterator is Traits::Curve_2, and the
 *      value-type of OutputIterator is Traits::X_monotone_curve_2.
 */
template <class Traits, class CurveInputIterator, class OutputIterator>
OutputIterator get_subcurves (CurveInputIterator curves_begin,
			      CurveInputIterator curves_end,
			      OutputIterator subcurves,
            Traits &tr,
			      bool mult_overlaps = false)
{
  // Define the sweep-line types:
  typedef Sweep_line_subcurves_visitor<Traits, OutputIterator>  Visitor;
  typedef Sweep_line_2<Traits,
                       Sweep_line_subcurves_visitor<Traits,
                                                    OutputIterator> >
                                                                Sweep_line;

  // Perform the sweep and obtain the subcurves.
  Visitor     visitor (subcurves, mult_overlaps, &tr);
  Sweep_line  sweep_line (&tr, &visitor);
  visitor.sweep(curves_begin, curves_end);

  return (visitor.get_output_iterator());
}

/*!
 * Determine if there occurs an intersection between any pair of curves in
 * a given range.
 * \param begin An input iterator for the first curve in the range.
 * \param end A input past-the-end iterator for the range.
 * \return (true) if any pair of curves intersect; (false) otherwise.
 */
template <class Traits, class CurveInputIterator>
bool do_curves_intersect (CurveInputIterator curves_begin,
                          CurveInputIterator curves_end,
                          Traits &tr)
{
  // Define the sweep-line types:
  typedef Sweep_line_do_curves_x_visitor<Traits>      Visitor;
  typedef Sweep_line_2<Traits,
                       Sweep_line_do_curves_x_visitor<Traits> >     
                                                      Sweep_line ;
  
  // Perform the sweep and obtain the subcurves.
  Visitor     visitor(&tr);
  Sweep_line  sweep_line (&tr, &visitor);
  visitor.sweep(curves_begin, curves_end);
  
  return (visitor.found_x());
}

CGAL_END_NAMESPACE

#endif
