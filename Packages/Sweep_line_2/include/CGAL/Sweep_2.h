// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.3-I-44 $
// release_date  : $CGAL_Date: 2001/03/09 $
//
// file          : include/CGAL/Sweep_2.h
// package       : arr (1.87)
// maintainer    : Tali Zvi <talizvi@post.tau.ac.il>
// source        : 
// revision      : 
// revision_date : 
// author(s)     : Tali Zvi <talizvi@post.tau.ac.il>
//
//
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// Chapter       : 
// ======================================================================
#ifndef CGAL_SWEEP_2_H
#define CGAL_SWEEP_2_H

#include <string>
#include "Sweep_line_tight_2.h"
#include <CGAL/Sweep_line_2/Sweep_line_event.h>
#include <CGAL/Sweep_line_2/Sweep_line_subcurve.h>

CGAL_BEGIN_NAMESPACE

/*!
  @class  Sweep_2

  Sweep_2 is a wrapper class to Sweep_line_tight_2. It is used as an interface
  to the Sweep algorithm.

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
  the input curves efficiently. \sa Pmwx_aggregate_insert_tight.

*/

template <class CurveInputIterator,  class SweepLineTraits_2>
class Sweep_2 : public Sweep_line_tight_2<
   CurveInputIterator, 
   SweepLineTraits_2, 
   Sweep_line_event<SweepLineTraits_2, Sweep_line_subcurve<SweepLineTraits_2> >,
   Sweep_line_subcurve<SweepLineTraits_2> >
{

public:
  typedef SweepLineTraits_2 Traits;
  typedef Sweep_line_tight_2<
    CurveInputIterator, 
    SweepLineTraits_2, 
    Sweep_line_event<Traits, Sweep_line_subcurve<Traits> >,
    Sweep_line_subcurve<Traits> > Base;

  Sweep_2() : Base() {}
  Sweep_2(Traits *t) : Base(t) {}
  virtual ~Sweep_2() {}

  /*!
   *  Given a container of curves, this function returns a list of curves
   *  that are created by intersecting the input curves.
   *  \param curves_begin the input iterator that points to the first curve 
   *                      in the range.
   *  \param curves_end the input past-the-end iterator of the range.
   *  \param subcurves an iterator to the first curve in the range
   *                   of curves created by intersecting the input curves.
   *  \param overlapping indicates whether overlapping curves should be 
   *                   reported once or multiple times. If false, the 
   *                   overlapping curves are reported once only.
   */
  template <class OutpoutIterator>
  void  get_subcurves(CurveInputIterator begin, CurveInputIterator end, 
		      OutpoutIterator subcurves, bool overlapping = false)
  { 
    Base::get_subcurves(begin, end, subcurves, overlapping);
  }

  /*!
   *  Given a range of curves, this function returns a list of points 
   *  that are the intersection points of the curves.
   *  The intersections are calculated using the sweep algorithm.
   *  \param curves_begin the input iterator that points to the first curve 
   *                      in the range.
   *  \param curves_end the input past-the-end iterator of the range.
   *  \param subcurves an iterator to the first curve in the range
   *                   of curves created by intersecting the input curves.
   *  \param endpoints if true, the end points of the curves are reported
   *                   as intersection points. Defaults to true.
   *  \param overlapping indicates whether there are overlapping curves
   *                     in the input range. Defaults to false.
   */
  template <class OutpoutIterator>
  void  get_intersection_points(CurveInputIterator begin, 
                                CurveInputIterator end, 
                                OutpoutIterator points,
                                bool includeEndPoints = true)
  { 
    Base::get_intersection_points(begin, end, points, includeEndPoints);
  }

 /*!
  *  Given a range of curves, this function returns an iterator 
  *  to the beginning of a range that contains the list of curves 
  *  for each intersection point between any two curves in the 
  *  specified range.
  *  The intersections are calculated using the sweep algorithm.
  *  \param curves_begin the input iterator that points to the first curve 
  *                      in the range.
  *  \param curves_end the input past-the-end iterator of the range.
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
    Base::get_intersecting_curves(begin, end, intersecting_curves, endpoints);
  }

};

CGAL_END_NAMESPACE

#endif // CGAL_SWEEP_2_H
