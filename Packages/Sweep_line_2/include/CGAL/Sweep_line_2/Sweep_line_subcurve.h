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
// file          : include/CGAL/Sweep_line_event.h
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

template<class SweepLineTraits_2>
class Sweep_line_subcurve
{
public:
  typedef SweepLineTraits_2 Traits;
  typedef typename Traits::Point_2 Point_2;
  typedef typename Traits::Curve_2 Curve_2;
  typedef typename Traits::X_curve_2 X_curve_2;

  Sweep_line_subcurve(X_curve_2 &curve, Point_2 *reference, 
		      SweepLineTraits_2 *traits);

  const X_curve_2 &getCurve() const { 
    return m_curve;
  }

  void setCurve(X_curve_2 &curve) {
    m_curve = curve;
  }

  const Point_2 *getReferencePoint() const { 
    return m_referencePoint;
  }

  const Point_2 &getLastPoint()  const { return m_lastPoint; }
  void setLastPoint(const Point_2 &point) { m_lastPoint = point; }

  const X_curve_2 &getLastCurve() const { return m_lastCurve; }
  void setLastCurve(const X_curve_2 &cv) { m_lastCurve = cv; }

  bool isSourceLeftToTarget() const { return m_isRightSide; }
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

#ifndef NDEBUG
  void Print() const;
#endif

private:
    
  Traits *m_traits;
  X_curve_2 m_curve;
  Point_2 *m_referencePoint;
  Point_2 m_lastPoint;
  X_curve_2 m_lastCurve;
  bool m_isRightSide;
  Point_2 m_source;
  Point_2 m_target;

public:
#ifndef NDEBUG
  int id;
#endif
};

template<class SweepLineTraits_2>
inline Sweep_line_subcurve<SweepLineTraits_2>::
Sweep_line_subcurve(X_curve_2 &curve, Point_2 *reference, 
		    SweepLineTraits_2 *traits)  : m_traits(traits)
{
  m_curve = curve;
  m_referencePoint = reference;
  m_source = traits->curve_source(curve);
  m_target = traits->curve_target(curve);
  if ( traits->compare_x(m_source, m_target)  == LARGER )
  {
    m_lastPoint =  m_target; //traits->curve_target(curve);
    m_isRightSide = false;
  }
  else
  {
    m_lastPoint =  m_source; //traits->curve_source(curve);
    m_isRightSide = true;
  }
}
#ifndef NDEBUG
template<class SweepLineTraits_2>
void 
Sweep_line_subcurve<SweepLineTraits_2>::
Print() const
{
  std::cout << "Curve " << id << "  (" << m_curve << ") "
	    << "last P = (" << m_lastPoint << ")" << std::endl;
  
}
#endif
CGAL_END_NAMESPACE

#endif // CGAL_SWEEP_LINE_SUBCURVE_H

