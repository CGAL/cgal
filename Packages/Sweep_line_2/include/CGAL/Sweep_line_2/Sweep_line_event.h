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
#ifndef CGAL_SWEEP_LINE_EVENT_H
#define CGAL_SWEEP_LINE_EVENT_H

#include <CGAL/Sweep_line_2/Sweep_line_subcurve.h>
#include <CGAL/Sweep_line_2/Sweep_line_functors.h>
#include <vector>
#include <set>

CGAL_BEGIN_NAMESPACE

template<class SweepLineTraits_2>
class Sweep_line_event
{
public:
  typedef SweepLineTraits_2 Traits;
  typedef typename Traits::X_curve_2 X_curve_2;
  typedef typename Traits::Point_2 Point_2;

  typedef Sweep_line_subcurve<Traits> SubCurve;
  typedef typename std::list<SubCurve *> SubcurveContainer;
  typedef typename SubcurveContainer::iterator SubCurveIter;

  typedef Status_line_curve_less_functor<Traits> StatusLineCurveLess;
  typedef std::set<SubCurve*, StatusLineCurveLess> StatusLine;
  typedef typename StatusLine::iterator StatusLineIter;

  /*! Constructor */
  Sweep_line_event(const Point_2 &point, Traits *traits) { 
    m_leftCurves = new SubcurveContainer();
    m_rightCurves = new SubcurveContainer();
    m_point = point;
    m_traits = traits;
    m_isInitialized = false;
    m_verticalCurve = 0;  
    m_isInternalIntersectionPoint = false;
  }

  /*! Destructor. Deletes the lists of curves, without deleting the 
      curves themselves. 
  */
  ~Sweep_line_event() {
    delete m_leftCurves;
    delete m_rightCurves;
  }


  /*! Adds a new curve to the event. The curve is added only to the list/s
      in which it is defined (left or/and right).
      The event point has to be either the source or the target of the curve.
      @param curve  a pointer to the curve.
  */
  void addCurve(SubCurve *scurve)
  {
    const X_curve_2 &curve = scurve->getCurve();
    const Point_2 &source = m_traits->curve_source(curve);
    const Point_2 &target = m_traits->curve_target(curve);
    
    if ( m_traits->curve_is_vertical(curve) ) 
    {
      m_verticalCurve = scurve;

    } else 
    {
      const Point_2 *rel = &(source);
      if ( m_traits->point_is_same(m_point, source) )
	rel = &(target);
      
      if ( m_traits->compare_x(m_point, *rel) == LARGER ) {
	addCurveToLeft1(scurve);
      } else {
	addCurveToRight(scurve);
      }
    }
  }

  /*! Adds a new curve that is defined to the right of the event point.
      The insertion is performed so that the curves remain srted by their Y 
      values to theright of the event.
      If the curve is already in the list of curves, it is removed and 
      re-inserted. This way the curves remain sorted.
      @param curve  a pointer to the curve.
  */
  void addCurveToLeft(SubCurve *curve, const Point_2 &ref) 
  {
    if (m_leftCurves->empty())
      m_leftCurves->push_back(curve);
    else 
    {
      SubCurveIter iter = m_leftCurves->begin();
      while ( iter != m_leftCurves->end() ) {
	if ( m_traits->curve_is_same((*iter)->getCurve(), curve->getCurve())) {
	  m_leftCurves->erase(iter);
	  break;
	}
	++iter;
      }
      iter = m_leftCurves->begin();
      while ( iter != m_leftCurves->end() &&
	      m_traits->curve_compare_at_x_right(curve->getCurve(),
						 (*iter)->getCurve(), 
						 ref)
	      == LARGER)
      {
	++iter;
      }

      if ( iter == m_leftCurves->end() ||
	   !m_traits->curve_is_same((*iter)->getCurve(), curve->getCurve()))
      {
	m_leftCurves->insert(iter, curve);
      } 
    }
    if ( !curve->isEndPoint(m_point) )
      m_isInternalIntersectionPoint = true;
  }

  /*! Adds a new curve that is defined to the right of the event point.
      The insertion is performed so that the curves remain srted by their Y 
      values to theright of the event.
      If the curve is already in the list of curves, it is removed and 
      re-inserted. This way the curves remain sorted.
      @param curve  a pointer to the curve.
  */
  void addCurveToLeft1(SubCurve *curve) 
  {
    if ( !m_isInitialized ) 
    {
      if ( curve->isSourceLeftToTarget())
	m_rightmostPointToLeft = curve->getSource();
      else
	m_rightmostPointToLeft = curve->getTarget();
      m_isInitialized = true;
    }

    UpdateRightmostPoint(curve);
    if (m_leftCurves->empty())
      m_leftCurves->push_back(curve);
    else 
    {
      SubCurveIter iter = m_leftCurves->begin();
      while ( iter != m_leftCurves->end() ) {
	if ( m_traits->curve_is_same((*iter)->getCurve(), curve->getCurve())) {
	  m_leftCurves->erase(iter);
	  break;
	}
	++iter;
      }
      iter = m_leftCurves->begin();
      while ( iter != m_leftCurves->end() &&
	      m_traits->curve_compare_at_x_right(curve->getCurve(),
						 (*iter)->getCurve(), 
						 m_rightmostPointToLeft)
	      == LARGER)
      {
	++iter;
      }

      if ( iter == m_leftCurves->end() ||
	   !m_traits->curve_is_same((*iter)->getCurve(), curve->getCurve()))
      {
	m_leftCurves->insert(iter, curve);
      }
    }
  }

  /*! Adds a new curve that is defined to the right of the event point.
      The insertion is performed so that the curves remain srted by their Y 
      values to theright of the event.
      @param curve  a pointer to the curve.
  */
  void addCurveToRight(SubCurve *curve) 
  {
    if (m_rightCurves->empty())
      m_rightCurves->push_back(curve);
    else 
    {
      SubCurveIter iter = m_rightCurves->begin();
      while ( iter != m_rightCurves->end() &&
	      m_traits->curve_compare_at_x_right(curve->getCurve(),
						 (*iter)->getCurve(), 
						 m_point)
	      == LARGER)
      {
	++iter;
      }
      if ( iter == m_rightCurves->end() ||
	   !m_traits->curve_is_same((*iter)->getCurve(), curve->getCurve()))
      {
	m_rightCurves->insert(iter, curve);
      }
    }
    if ( !curve->isEndPoint(m_point) )
      m_isInternalIntersectionPoint = true;
  }

  /*! Returns an iterator to the first curve to the left of the event */
  SubCurveIter leftCurvesBegin() {
    return m_leftCurves->begin();
  }

  /*! Returns an iterator to the one past the last curve to the left 
      of the event */
  SubCurveIter leftCurvesEnd() {
    return m_leftCurves->end();
  }

  /*! Returns an iterator to the first curve to the right of the event */
  SubCurveIter rightCurvesBegin() {
    return m_rightCurves->begin();
  }

  /*! Returns an iterator to the one past the last curve to the right 
      of the event */
   SubCurveIter rightCurvesEnd() {
    return m_rightCurves->end();
  }

  /*! Returns the number of intersecting curves that are defined
      to the right of the event point. */
  int getNumRightCurves() {
    return m_rightCurves->size();
  }

  /*! Returns true if at least one intersecting curve is defined to 
      the left of the point. */
  bool hasLeftCurves() {
    return !m_leftCurves->empty();
  }

  /*! Returns the actual point of the event */
  const Point_2 &getPoint() {
    return m_point;
  }

  bool doesContainVerticalCurve() {
    return m_verticalCurve != 0;
  }
  SubCurve *getVerticalCurve() {
    return m_verticalCurve;
  }

  void addVerticalCurveXPoint(const Point_2 &p) {
    if ( m_verticalCurveXPoints.empty() ) 
    {
      m_verticalCurveXPoints.push_back(p); 
      return;
    }
    if (!m_traits->point_is_same(p, m_verticalCurveXPoints.back())) {
      m_verticalCurveXPoints.push_back(p);
    }
  }

  std::list<Point_2> &getVerticalCurveXPointList() {
    return m_verticalCurveXPoints;
  }

  bool isInternalIntersectionPoint() {
    return m_isInternalIntersectionPoint;
  }
  void markInternalIntersectionPoint() {
    m_isInternalIntersectionPoint = true;
  }

#ifndef NDEBUG
  void Print();
#endif
 
private:

  void UpdateRightmostPoint(SubCurve *curve)
  {
    if ( curve->isSourceLeftToTarget())
    {
      if ( curve->isTarget(m_point) )
	if ( m_traits->compare_x(m_point, m_rightmostPointToLeft) == LARGER )
	  m_rightmostPointToLeft = m_point;
    } else
    {
      if ( curve->isSource(m_point) )
	if ( m_traits->compare_x(m_point, m_rightmostPointToLeft) == LARGER )
	  m_rightmostPointToLeft = m_point;
    }
  }

  /*! A pointer to a traits class */
  Traits *m_traits;

  /*! A list of curves on the left side of the event, sorted by their y value
      to the left of the point */
  SubcurveContainer *m_leftCurves;

  /*! A list of curves on the right side of the event, sorted by their y value
      to the right of the point */
  SubcurveContainer *m_rightCurves;

  /*! The point of the event */
  Point_2 m_point;

  /*! The rightmost curve-end point that is to the left of the event
      point. This point is used as a reference point when curves are compared
      only when new curves are added at the initialization stage.
  */
  Point_2 m_rightmostPointToLeft;

  /*! An indication whether this event has been initialized. The event is
      initialized after the first curve has been added to the left of the 
      event. 
  */
  bool m_isInitialized;

  /*! a pointer to a vertical curve going through this event */
  SubCurve *m_verticalCurve;

  /*! a list of intersection points on the vertical curve */
  std::list<Point_2> m_verticalCurveXPoints;

  /*! a flag that inidcates whether the event is an "interior" intersection 
      point.
  */
  bool m_isInternalIntersectionPoint;

#ifndef NDEBUG
public:
  int id;
#endif
  
};

#ifndef NDEBUG
template<class SweepLineTraits_2>
void 
Sweep_line_event<SweepLineTraits_2>::
Print() 
{
  std::cout << "\tEvent id: " << id << "\n" ;
  std::cout << "\tLeft curves: \n" ;
  for ( SubCurveIter iter = m_leftCurves->begin() ;
	iter != m_leftCurves->end() ; ++iter )
  {
    std::cout << "\t";
    (*iter)->Print();
    std::cout << "\n";
  }
  std::cout << std::endl;
  std::cout << "\tRight curves: \n" ;
  for ( SubCurveIter iter = m_rightCurves->begin() ;
	iter != m_rightCurves->end() ; ++iter )
  {
    std::cout << "\t";
    (*iter)->Print();
    std::cout << "\n";
  }
  std::cout << std::endl;
}
#endif // NDEBUG

CGAL_END_NAMESPACE

#endif // CGAL_SWEEP_LINE_EVENT_H
