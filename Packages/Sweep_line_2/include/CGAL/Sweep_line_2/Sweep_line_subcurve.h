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
#ifndef CGAL_SWEEP_LINE_SUBCURVE_H
#define CGAL_SWEEP_LINE_SUBCURVE_H

#include <vector>
#include <set>
#include <CGAL/Sweep_line_2/Sweep_line_functors.h>
#include <CGAL/Sweep_line_2/Sweep_line_event.h>
#include <CGAL/assertions.h>

CGAL_BEGIN_NAMESPACE

/*! @class Sweep_line_subcurve
 *
 * This is a wrapper class to Curve_2 in the traits class, that contains
 * data that is used when applying the sweep algorithm on a set of curves.
 *
 * The information contained in this class is:
 * - the curve itself
 * - two points which are the source and target of the curve. We keep 
 *   the points in order to avoid many calls to the source() and 
 *   target() methods of the traits class 
 * - an indication for the direction of the curve (source point 
 *   is left or right to the target point). 
 * - a reference point that is used when comparing the y values of 
 *   any two curves. Since the curves are inserted in to a balanced 
 *   tree, and at any given time they are sorted on the status line, 
 *   and since their order may change, depending on the position of 
 *   the status line, we need to be able to compare the curves 
 *   relative to a point that will produce a correct answer.
 * - a reference to the last event point on the curve that was already 
 *   handled and also the curve that is the portion of the original 
 *   curve that is to the right of the last event point mentioned. 
 *   This is stored to avoid unneccesary splits of the curve.
 *
 */

template<class SweepLineTraits_2>
class Sweep_line_subcurve
{
public:
  typedef SweepLineTraits_2 Traits;
  typedef typename Traits::Point_2 Point_2;
  typedef typename Traits::Curve_2 Curve_2;

  typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;

  typedef Sweep_line_subcurve<Traits> Self;
  typedef Status_line_curve_less_functor<Traits, Self> StatusLineCurveLess;

  typedef std::set<Self*, StatusLineCurveLess> StatusLine;
  typedef typename StatusLine::iterator StatusLineIter;

  typedef Sweep_line_event<Traits, Self> Event;

  Sweep_line_subcurve(){}

  Sweep_line_subcurve(X_monotone_curve_2 &curve ,SweepLineTraits_2 *traits);

  void init(X_monotone_curve_2 &curve ,SweepLineTraits_2 *traits)
  {
    m_traits = traits;
    m_curve = curve;
    m_source = traits->curve_source(curve);
    m_target = traits->curve_target(curve);
    Comparison_result res = traits->compare_xy(m_source, m_target);
    if ( res  == LARGER )
    {
      m_lastPoint = m_target; 
      m_isRightSide = false;
    }
    else //if ( res == SMALLER )
    { 
      CGAL_assertion(res == SMALLER); //curves cannot be a degenerate point
      m_lastPoint = m_source; 
      m_isRightSide = true;
    }
    m_lastCurve = curve;
  }

  virtual ~Sweep_line_subcurve() {}

  /*int getId() const {
    return m_id;
  }*/

  /*!
    @return a reference to the curve 
  */
  const X_monotone_curve_2 &get_curve() const { 
    return m_curve;
  }

  ///*! @return  the pointer to the reference point */
  //const Point_2 *get_reference_point() const { 
  //  return m_referencePoint;
  //}

  /*! 
    @return a reference to the rightmost intersection point 
  */
  const Point_2 &get_last_point()  const { 
    return m_lastPoint; 
  }

  /*! 
    Updates the rightmost intersection point.
    @param point a reference to the point
   */
  void set_last_point(const Point_2 &point) { 
    m_lastPoint = point; 
  }

  /*!
    @return a const reference to the last intersecing curve so far
  */
  const X_monotone_curve_2 &get_last_curve() const { 
    return m_lastCurve; 
  }

  /*!
    @return a reference to the last intersecing curve so far
  */
  X_monotone_curve_2 &get_last_curve_ref() { 
    return m_lastCurve; 
  }

  /*! 
    updates the last intersecting curve so far.
    @param cv a reference to the curve
  */
  void set_last_curve(const X_monotone_curve_2 &cv) { 
    m_lastCurve = cv; 
  }

  const X_monotone_curve_2 &get_last_subcurve() const { 
    return m_lastSubCurve; 
  }

 /* X_monotone_curve_2 &get_last_subcurve_ref()  { 
    return m_lastSubCurve; 
  }*/

  void set_last_subcurve(const X_monotone_curve_2 &cv) { 
    m_lastSubCurve = cv; 
  }


  bool is_source_left_to_target() const { 
    return m_isRightSide; 
  }

  bool is_source(const Point_2 &p) { 
    return m_traits->point_equal(p, m_source);
  }

  bool is_target(const Point_2 &p) { 
    return m_traits->point_equal(p, m_target);
  }

  /*! returns true if the specified point is the source or the target
      of the curve. Returns false otherwise.
  */
  bool is_end_point(const Point_2 &p) { 
    return is_target(p) || is_source(p);
  }

  /*! returns true if the last point is an end point and the specified
      point is an end point. Otherwise returns false;
  */
  bool is_unsplit_curve(const Point_2 &p) {
    if ( is_end_point(p) && is_end_point(m_lastPoint) )
      return true;
    return false;
  }

  const Point_2 &get_source() const {
    return m_source;
  }

  const Point_2 &get_target() const {
    return m_target;
  }

  bool is_left_end(const Point_2 &p) {
    if ( is_source_left_to_target() && is_source(p) )
      return true;
    if ( !is_source_left_to_target() && is_target(p) )
      return true;
    return false;
  }

  bool is_right_end(const Point_2 &p) {
    if ( is_source_left_to_target() && is_target(p) )
      return true;
    if ( !is_source_left_to_target() && is_source(p) )
      return true;
    return false;
  }

  /*bool is_bottom_end(const Point_2 &p) {
    CGAL_assertion(m_traits->curve_is_vertical(m_curve)==true);
    return is_left_end(p);
  }

  bool is_top_end(const Point_2 &p) {
    CGAL_assertion(m_traits->curve_is_vertical(m_curve)==true);
    return is_right_end(p);
  }*/

  const Point_2 &get_right_end() const {
    if ( is_source_left_to_target() )
      return m_target;
    return m_source;
  }

  const Point_2 &get_left_end() const {
    if ( is_source_left_to_target() )
      return m_source;
    return m_target;
  }

  /*const Point_2 &get_top_end() const {
    CGAL_assertion(m_traits->curve_is_vertical(m_curve)==true);
    return get_right_end();
  }

  const Point_2 &get_bottom_end() const {
    CGAL_assertion(m_traits->curve_is_vertical(m_curve)==true);
    return get_left_end();
  }*/

  // returns true if the point is in the range of the curve and is not
  // one of its ends
  bool is_point_in_range(const Point_2 &p)
  {
    if (! m_traits->point_in_x_range(m_curve, p) ||
        m_traits->curve_compare_y_at_x(p, m_curve) != EQUAL)
      return false;
    if ( is_end_point(p) )
      return false;
    return true;
  }


  template <class StatusLineIter>
  void set_hint(StatusLineIter hint) 
  {
    m_hint = hint;
  }

  StatusLineIter get_hint() const 
  {
    return m_hint;
  }

#ifndef NDEBUG
  void Print() const;
#endif

private:

  //int m_id;

  /*! a pointer to the traits object */
  Traits *m_traits;

  /*! thecurve */
  X_monotone_curve_2 m_curve;

  /* a pointer to a point that is used as a reference point when two 
     curves are compared. This is used when inserting and erasing 
     curves from the status line. */
 // Point_2 *m_referencePoint;

  /*! the rightmost point handled so far on the curve. It is initialized 
    to the left end of the curve and is updated with every intersection 
    point on the curve. */
  Point_2 m_lastPoint;

  /*! the portion of the curve to the right of the last event point 
      on the curve */
  X_monotone_curve_2 m_lastCurve;

  /*! the last subcurve that was reported */
  X_monotone_curve_2 m_lastSubCurve;

  /*! true if the source of the curve is to the left of the target. */
  bool m_isRightSide;

  /*! the source of the curve */
  Point_2 m_source;

  /*! the target of the curve */
  Point_2 m_target;

  /*! */
  StatusLineIter m_hint;

};

template<class SweepLineTraits_2>
inline Sweep_line_subcurve<SweepLineTraits_2>::
Sweep_line_subcurve( X_monotone_curve_2 &curve,
                    SweepLineTraits_2 *traits)  :  m_traits(traits)
{
  m_curve = curve;
 // m_referencePoint = reference;
  m_source = traits->curve_source(curve);
  m_target = traits->curve_target(curve);
  Comparison_result res = traits->compare_xy(m_source, m_target);
  if ( res  == LARGER )
  {
    m_lastPoint = m_target; 
    m_isRightSide = false;
  }
  else //if ( res == SMALLER )
  { 
    CGAL_assertion(res == SMALLER); //curves cannot be a degenerate point
    m_lastPoint = m_source; 
    m_isRightSide = true;
  }
  m_lastCurve = curve;
}

#ifndef NDEBUG
template<class SweepLineTraits_2>
void 
Sweep_line_subcurve<SweepLineTraits_2>::
Print() const
{
  std::cout << "Curve " << this << "  (" << m_curve << ") "
            << "last P = (" << m_lastPoint << ")" << std::endl;
  
}

#endif

CGAL_END_NAMESPACE

#endif
