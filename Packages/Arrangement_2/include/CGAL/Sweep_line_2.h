// Copyright (c) 1997  Tel-Aviv University (Israel).
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
// Author(s)     : Tali Zvi <talizvi@post.tau.ac.il>
#ifndef CGAL_SWEEP_2_H
#define CGAL_SWEEP_2_H

#include <CGAL/Sweep_line_2/Sweep_line_2_impl.h>
#include <CGAL/Sweep_line_2/Sweep_line_event.h>
#include <CGAL/Sweep_line_2/Sweep_line_subcurve.h>
#include <CGAL/Sweep_line_2/Sweep_line_points_visitor.h>
#include <CGAL/Sweep_line_2/Sweep_line_subcurves_visitor.h>
#include <CGAL/Sweep_line_2/Sweep_line_do_curves_x_visitor.h>





CGAL_BEGIN_NAMESPACE


 /*!
  *  Given a range of curves, this function returns a list of points 
  *  that are the intersection points of the curves.
  *  The intersections are calculated using the sweep algorithm.
  *  \param begin the input iterator that points to the first curve 
  *                   in the range.
  *  \param end the input past-the-end iterator of the range.
  *  \param points an iterator to the first point in the range
  *                   of points created by intersecting the input curves.
  *  \param includeEndPoints if true, the end points of the curves are 
  *                   reported as intersection points. Defaults to true.
  */
template <class Traits, class CurveInputIterator, class OutputIterator>
OutputIterator get_intersection_points(CurveInputIterator curves_begin,
                                       CurveInputIterator curves_end,
                                       OutputIterator points,
                                       bool endpoints = true,
                                       Traits traits  = Traits() )
{
  typedef Sweep_line_points_visitor<Traits,OutputIterator> Visitor;
  Visitor visitor(points, endpoints);
  
  typedef Sweep_line_subcurve<Traits>                         Subcurve;
  
  typedef Sweep_line_event<Traits, Subcurve>            Event;
  typedef Sweep_line_2_impl< Traits,
                             Event,
                             Subcurve,
                             Sweep_line_points_visitor<Traits,OutputIterator>,
                             CGAL_ALLOCATOR(int) > Sweep_line ;
  Sweep_line sweep_object(&traits, &visitor);
  sweep_object.init(curves_begin, curves_end);
  sweep_object.sweep();

  return visitor.get_output_iterator();
}





 /*!
  *  Given a container of curves, this function returns a list of curves
  *  that are created by intersecting the input curves.
  *  \param begin the input iterator that points to the first curve 
  *                      in the range.
  *  \param end the input past-the-end iterator of the range.
  *  \param subcurves an iterator to the first curve in the range
  *                   of curves created by intersecting the input curves.
  *  \param overlapping indicates whether overlapping curves should be 
  *                   reported once or multiple times. If false, the 
  *                   overlapping curves are reported once only.
  */
template <class Traits, class CurveInputIterator, class OutputIterator>
OutputIterator get_subcurves(CurveInputIterator curves_begin,
                             CurveInputIterator curves_end,
                             OutputIterator subcurves,
                             bool overlapping = false,
                             Traits traits = Traits())
{
  typedef Sweep_line_subcurves_visitor<Traits, OutputIterator>    Visitor;
  Visitor visitor(subcurves, overlapping);


  typedef Sweep_line_subcurve<Traits>                         Subcurve;
  
  typedef Sweep_line_event<Traits, Subcurve>            Event;
  typedef Sweep_line_2_impl< Traits,
                             Event,
                             Subcurve,
                             Sweep_line_subcurves_visitor<Traits,OutputIterator>,
                             CGAL_ALLOCATOR(int) > Sweep_line ;

  Sweep_line sweep_object(&traits, &visitor);
  sweep_object.init(curves_begin, curves_end);
  sweep_object.sweep();
  return visitor.get_output_iterator();
}




                                                                  



/*!
  *  Given a range of curves, this function returns an iterator 
  *  to the beginning of a range that contains the list of curves 
  *  for each intersection point between any two curves in the 
  *  specified range.
  *  The intersections are calculated using the sweep algorithm.
  *  \param begin the input iterator that points to the first curve 
  *                      in the range.
  *  \param end the input past-the-end iterator of the range.
  *  \param intersecting_curves an iterator to the output
  *  \param endpoints if true, the end points of the curves are reported
  *                   as intersection points. Defaults to true.
  */
template <class Traits, class CurveInputIterator, class OutputIterator>
OutputIterator get_intersecting_curves(CurveInputIterator curves_begin,
                                       CurveInputIterator curves_end,
                                       OutputIterator intersecting_curves,
                                       bool endpoints = true,
                                       Traits traits = Traits())
{

  //TODO: implement !!!
  return intersecting_curves;
}


template <class Traits, class CurveInputIterator>
bool do_curves_intersect(CurveInputIterator curves_begin,
                         CurveInputIterator curves_end,
                         Traits traits = Traits())
{
  typedef Sweep_line_do_curves_x_visitor<Traits>  Visitor;
  Visitor visitor;

  typedef Sweep_line_subcurve<Traits>                         Subcurve;
  
  typedef Sweep_line_event<Traits, Subcurve>            Event;
  typedef Sweep_line_2_impl< Traits,
                             Event,
                             Subcurve,
                             Sweep_line_do_curves_x_visitor<Traits>,
                             CGAL_ALLOCATOR(int) > Sweep_line ;
  
  Sweep_line sweep_object(&traits, &visitor);
  sweep_object.init(curves_begin, curves_end);
  sweep_object.sweep();
  return visitor.found_x();
}





/*!
  @class  Sweep_line_2

  Sweep_line_2 is a wrapper class to Sweep_line_2_impl. It is used as an 
  interface to the Sweep algorithm.

  The sweep suupport the following degenerate cases:
  - non x-monotone curves
  - vertical segments
  - multiple (more then two) segments intersecting at one point
  - curves beginning and ending on other curves.
  - overlapping curves

  There are two main functionalities supported by this algorithm:
  1. calculate the non intersecting curves that are a product of 
     intersecting a set of input curves
  2. calculate all the intersection points between the curves specified

  An extension of this algorithm is used to produce a planar map containing
  the input curves efficiently. 
  \sa Pmwx_aggregate_insert.

*/


template <class CurveInputIterator,  class SweepLineTraits_2>
class Sweep_line_2 
{

public:

  typedef SweepLineTraits_2 Traits;
 
  Sweep_line_2() {}

  Sweep_line_2(Traits *t) {}

  virtual ~Sweep_line_2() {}
 

  /*!
   *  Given a container of curves, this function returns a list of curves
   *  that are created by intersecting the input curves.
   *  \param begin the input iterator that points to the first curve 
   *                      in the range.
   *  \param end the input past-the-end iterator of the range.
   *  \param subcurves an iterator to the first curve in the range
   *                   of curves created by intersecting the input curves.
   *  \param overlapping indicates whether overlapping curves should be 
   *                   reported once or multiple times. If false, the 
   *                   overlapping curves are reported once only.
   */
  template <class OutputIterator>
  void  get_subcurves(CurveInputIterator begin, CurveInputIterator end, 
		      OutputIterator subcurves, bool overlapping = false)
  { 
    CGAL::get_subcurves(begin, end, subcurves, overlapping, Traits());
  }

  /*!
   *  Given a range of curves, this function returns a list of points 
   *  that are the intersection points of the curves.
   *  The intersections are calculated using the sweep algorithm.
   *  \param begin the input iterator that points to the first curve 
   *                   in the range.
   *  \param end the input past-the-end iterator of the range.
   *  \param points an iterator to the first point in the range
   *                   of points created by intersecting the input curves.
   *  \param includeEndPoints if true, the end points of the curves are 
   *                   reported as intersection points. Defaults to true.
   */
  template <class OutputIterator>
  void  get_intersection_points(CurveInputIterator begin, 
                                CurveInputIterator end, 
                                OutputIterator points,
                                bool includeEndPoints = true)
  { 
    CGAL::get_intersection_points(begin,end,points,includeEndPoints, Traits());
  }

 /*!
  *  Given a range of curves, this function returns an iterator 
  *  to the beginning of a range that contains the list of curves 
  *  for each intersection point between any two curves in the 
  *  specified range.
  *  The intersections are calculated using the sweep algorithm.
  *  \param begin the input iterator that points to the first curve 
  *                      in the range.
  *  \param end the input past-the-end iterator of the range.
  *  \param intersecting_curves an iterator to the output
  *  \param endpoints if true, the end points of the curves are reported
  *                   as intersection points. Defaults to true.
  */
  template <class OutputIterator>
  void  get_intersecting_curves(CurveInputIterator begin, 
				CurveInputIterator end, 
				OutputIterator intersecting_curves,
				bool endpoints = true)
  { 
    CGAL::get_intersecting_curves(begin, end, intersecting_curves, endpoints,Traits());
  }

  bool do_curves_intersect(CurveInputIterator begin, CurveInputIterator end)
  {
    return CGAL::do_curves_intersect(begin, end, Traits());
  }

};

CGAL_END_NAMESPACE

#endif
