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
// file          : include/CGAL/Sweep_line_subcurve.h
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
#ifndef CGAL_SWEEP_LINE_SUBCURVE_H
#define CGAL_SWEEP_LINE_SUBCURVE_H

#include <vector>

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
  typedef typename Traits::X_curve_2 X_curve_2;

  Sweep_line_subcurve(int id, X_curve_2 &curve, Point_2 *reference, 
		      SweepLineTraits_2 *traits);


  int getId() const {
    return m_id;
  }

  /*!
    @return a reference to the curve 
  */
  const X_curve_2 &getCurve() const { 
    return m_curve;
  }

  /*! @return  the pointer to the reference point */
  const Point_2 *getReferencePoint() const { 
    return m_referencePoint;
  }

  /*! 
    @return a reference to the rightmost intersection point 
  */
  const Point_2 &getLastPoint()  const { 
    return m_lastPoint; 
  }

  /*! 
    Updates the rightmost intersection point.
    @param point a reference to the point
   */
  void setLastPoint(const Point_2 &point) { 
    m_lastPoint = point; 
  }

  /*!
    @return a reference to the last intersecing curve so far
  */
  const X_curve_2 &getLastCurve() const { 
    return m_lastCurve; 
  }
  /*! 
    updates the last intersecting curve so far.
    @param cv a reference to the curve
  */
  void setLastCurve(const X_curve_2 &cv) { 
    m_lastCurve = cv; 
  }

  bool isSourceLeftToTarget() const { 
    return m_isRightSide; 
  }

  bool isSource(const Point_2 &p) { 
    return m_traits->point_is_same(p, m_source);
  }

  bool isTarget(const Point_2 &p) { 
    return m_traits->point_is_same(p, m_target);
  }

  /*! returns true if the specified point is the source or the target
      of the curve. Returns false otherwise.
  */
  bool isEndPoint(const Point_2 &p) { 
    return isTarget(p) || isSource(p);
  }

  /*! returns true if the last point is an end point and the specified
      point is an end point. Otherwise returns false;
  */
  bool isUnsplitCurve(const Point_2 &p) {
    if ( isEndPoint(p) && isEndPoint(m_lastPoint) )
      return true;
    return false;
  }

  const Point_2 &getSource() const {
    return m_source;
  }

  const Point_2 &getTarget() const {
    return m_target;
  }

  bool isLeftEnd(const Point_2 &p) {
    if ( isSourceLeftToTarget() && isSource(p) )
      return true;
    if ( !isSourceLeftToTarget() && isTarget(p) )
      return true;
    return false;
  }

  bool isRightEnd(const Point_2 &p) {
    if ( isSourceLeftToTarget() && isTarget(p) )
      return true;
    if ( !isSourceLeftToTarget() && isSource(p) )
      return true;
    return false;
  }

  bool isBottomEnd(const Point_2 &p) {
    assert(m_traits->curve_is_vertical(m_curve)==true);
    return isLeftEnd(p);
  }

  bool isTopEnd(const Point_2 &p) {
    assert(m_traits->curve_is_vertical(m_curve)==true);
    return isRightEnd(p);
  }

  const Point_2 &getRightEnd() {
    if ( isSourceLeftToTarget() )
      return m_target;
    return m_source;
  }

  const Point_2 &getLeftEnd() {
    if ( isSourceLeftToTarget() )
      return m_source;
    return m_target;
  }

  const Point_2 &getTopEnd() {
    assert(m_traits->curve_is_vertical(m_curve)==true);
    return getRightEnd();
  }

  const Point_2 &getBottomEnd() {
    assert(m_traits->curve_is_vertical(m_curve)==true);
    return getLeftEnd();
  }

  // returns true if the point is in the range of the curve and is not
  // one of its ends
  bool isPointInRange(const Point_2 &p)
  {
    if (! m_traits->curve_is_in_x_range(m_curve, p) ||
	m_traits->curve_get_point_status(m_curve, p) != EQUAL)
      return false;
    if ( isEndPoint(p) )
      return false;
    return true;
  }

#ifndef NDEBUG
  void Print() const;
#endif

private:

  int m_id;

  /*! a pointer to the traits object */
  Traits *m_traits;

  /*! thecurve */
  X_curve_2 m_curve;

  /* a pointer to a point that is used as a reference point when two 
     curves are compared. This is used when inserting and erasing 
     curves from the status line. */
  Point_2 *m_referencePoint;

  /*! the rightmost point handled so far on the curve. It is initialized 
    to the left end of the curve and is updated with every intersection 
    point on the curve. */
  Point_2 m_lastPoint;

  /*! the portion of the curve to the right of the last event point 
      on the curve */
  X_curve_2 m_lastCurve;

  /*! true if the source of the curve is to the left of the target. */
  bool m_isRightSide;

  /*! the source of the curve */
  Point_2 m_source;

  /*! the target of the curve */
  Point_2 m_target;

};

template<class SweepLineTraits_2>
inline Sweep_line_subcurve<SweepLineTraits_2>::
Sweep_line_subcurve(int id, X_curve_2 &curve, Point_2 *reference, 
		    SweepLineTraits_2 *traits)  : m_id(id), m_traits(traits)
{
  m_curve = curve;
  m_referencePoint = reference;
  m_source = traits->curve_source(curve);
  m_target = traits->curve_target(curve);
  Comparison_result res = traits->compare_xy(m_source, m_target);
  if ( res  == LARGER )
  {
    m_lastPoint = m_target; 
    m_isRightSide = false;
  }
  else if ( res == SMALLER )
  {
    m_lastPoint = m_source; 
    m_isRightSide = true;

  }
}

#ifndef NDEBUG
template<class SweepLineTraits_2>
void 
Sweep_line_subcurve<SweepLineTraits_2>::
Print() const
{
  std::cout << "Curve " << m_id << "  (" << m_curve << ") "
	    << "last P = (" << m_lastPoint << ")" << std::endl;
  
}
#endif
CGAL_END_NAMESPACE

#endif // CGAL_SWEEP_LINE_SUBCURVE_H

