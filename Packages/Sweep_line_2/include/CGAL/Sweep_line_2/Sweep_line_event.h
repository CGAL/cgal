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

#include <CGAL/Sweep_line_2/Sweep_line_functors.h>
#include <list>
#include <set>

CGAL_BEGIN_NAMESPACE

template<class SweepLineTraits_2, class CurveWrap>
class Sweep_line_event
{
public:
  typedef SweepLineTraits_2 Traits;
  typedef typename Traits::X_curve_2 X_curve_2;
  typedef typename Traits::Point_2 Point_2;

  //typedef Sweep_line_subcurve<Traits> SubCurve;
  typedef CurveWrap SubCurve;
  typedef typename std::list<SubCurve *> SubcurveContainer;
  typedef typename SubcurveContainer::iterator SubCurveIter;

  typedef Status_line_curve_less_functor<Traits, SubCurve> StatusLineCurveLess;
  typedef std::set<SubCurve*, StatusLineCurveLess> StatusLine;
  typedef typename StatusLine::iterator StatusLineIter;

  typedef std::list<Point_2> VerticalXPointList;
  typedef VerticalXPointList::iterator VerticalXPointListIter; 

  typedef std::list<SubCurve *> VerticalCurveList;
  typedef VerticalCurveList::iterator VerticalCurveListIter;

  /*! Constructor */
  Sweep_line_event(const Point_2 &point, Traits *traits) :
    m_point(point), m_traits(traits), m_isInitialized(false),
    m_isInternalIntersectionPoint(false), m_containsOverlap(false)
  { 
    m_leftCurves = new SubcurveContainer();
    m_rightCurves = new SubcurveContainer();
  }

  /*! Destructor. Deletes the lists of curves, without deleting the 
      curves themselves. 
  */
  virtual ~Sweep_line_event() {
    delete m_leftCurves;
    delete m_rightCurves;
  }


  /*! Adds a new curve to the event. The curve is added only to the list/s
      in which it is defined (left or/and right).
      If the curve is vertical, it is added to the list of vertical curves.
      Precondition: The event point has to be either the source or the 
      target of the curve.
      @param curve  a pointer to the curve.
  */
  void addCurve(SubCurve *scurve)
  {
    const X_curve_2 &curve = scurve->getCurve();
    const Point_2 &source = m_traits->curve_source(curve);
    const Point_2 &target = m_traits->curve_target(curve);
    
    if ( m_traits->curve_is_vertical(curve) ) 
    {
      m_verticalCurves.push_back(scurve);

    } else 
    {
      const Point_2 *rel = &(source);
      if ( m_traits->point_is_same(m_point, source) )
	rel = &(target);
      
      if ( m_traits->compare_x(m_point, *rel) == LARGER ) {
	addCurveToLeft(scurve, m_rightmostPointToLeft, true);
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
      @pram ref a reference point to perform the compare by
      @param isInitStage true when thie method is called at the 
      initialization stage (in which case some extra tests are performed).
  */
  void addCurveToLeft(SubCurve *curve, const Point_2 &ref, 
		      bool isInitStage=false) 
  {
    if ( isInitStage )
    {
      if ( !m_isInitialized ) 
      {
	if ( curve->isSourceLeftToTarget())
	  m_rightmostPointToLeft = curve->getSource();
	else
	  m_rightmostPointToLeft = curve->getTarget();
	m_isInitialized = true;

      } else {
	UpdateRightmostPoint(curve);
      }
    } else if ( !curve->isEndPoint(m_point) ) {

      m_isInternalIntersectionPoint = true;
    }

    // now insert the curve at the right place...

    if (m_leftCurves->empty()) {
      m_leftCurves->push_back(curve);
      return;
    }

    SubCurveIter iter = m_leftCurves->begin();
    const X_curve_2 &cv = curve->getCurve();
    
    // look for the curve, and if exists, erase it.
    while ( iter != m_leftCurves->end() ) {
      if ( (*iter)->getId() ==  curve->getId()) {
	m_leftCurves->erase(iter);
	break;
      }
      ++iter;
    }
    
    // insert the curve so that the list remains sorted...
    Comparison_result res = SMALLER;
    iter = m_leftCurves->begin();

    while ( iter != m_leftCurves->end() &&
	    (res = m_traits->curve_compare_at_x_right(cv, (*iter)->getCurve(), 
						      ref )) == LARGER)
      ++iter;
    
    while ( iter != m_leftCurves->end() &&
	    res == EQUAL &&
	    curve->getId() > (*iter)->getId() )
    {
      m_containsOverlap = true;
      ++iter;
      if ( iter == m_leftCurves->end())
	break;
      res = m_traits->curve_compare_at_x_right(cv, (*iter)->getCurve(), ref);
    }
    
    // insert the curve. If the curve is already in the list, it is not added
    if ( iter == m_leftCurves->end() ||
	 (*iter)->getId() !=  curve->getId())
    {
      m_leftCurves->insert(iter, curve);
    } 
  }


  /*! Adds a new curve that is defined to the right of the event point.
      The insertion is performed so that the curves remain sorted by their Y 
      values to the right of the event.
      @param curve  a pointer to the curve.
  */
  void addCurveToRight(SubCurve *curve) 
  {
    if ( !curve->isEndPoint(m_point) )
      m_isInternalIntersectionPoint = true;

    if (m_rightCurves->empty()) {
      m_rightCurves->push_back(curve);
      return;
    }


    SubCurveIter iter = m_rightCurves->begin();
    Comparison_result res;
    while ( (res = m_traits->curve_compare_at_x_right(curve->getCurve(),
						      (*iter)->getCurve(), 
						      m_point)) == LARGER)
    {
      ++iter;
      if ( iter == m_rightCurves->end()) {
	m_rightCurves->insert(iter, curve);
	return;
      }
    }
    
    while ( res == EQUAL && curve->getId() > (*iter)->getId() )
    {
      m_containsOverlap = true;
      ++iter;
      if ( iter == m_rightCurves->end() ) {
	m_rightCurves->insert(iter, curve);
	return;
      }
      res = m_traits->curve_compare_at_x_right(curve->getCurve(),
					       (*iter)->getCurve(), 
					       m_point);
    }
    
    // insert the curve only if it is not already in...
    if ( (*iter)->getId() !=  curve->getId()) {
      m_rightCurves->insert(iter, curve);
    }
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

  /*! 
    @return returns true if at least one of the curves passign 
    through the event is vertical.
  */
  bool doesContainVerticalCurve() const {
    return !m_verticalCurves.empty();
  }


  /*! eturns the list of vertical curves passing through the event point.
    @return a reference to the list of curves.
  */
  VerticalCurveList &getVerticalCurves() {
    return m_verticalCurves;
  }


  /*! Insert a new intersection point on any of the vertical curves.
      The list of points is sorted by their y values.
      If the requireSort flag is true, the appripriate place in the list 
      is searched for. If not, the point is assumed to have the largest y 
      value, and is inserted at the end of the list. 
      If the pioint already exists, the point is nott inserted again.
      @param p a reference to the point
      @param requireSort false if the point is to be added at the end
      of the list.
  */
  void addVerticalCurveXPoint(const Point_2 &p, bool requireSort=false) {
    if ( m_verticalCurveXPoints.empty() ) 
    {
      m_verticalCurveXPoints.push_back(p); 
      return;
    }

    if ( !requireSort ) 
    {
      if (!m_traits->point_is_same(p, m_verticalCurveXPoints.back())) {
	m_verticalCurveXPoints.push_back(p);
      }
    } else
    {
      VerticalXPointListIter iter = m_verticalCurveXPoints.begin();
      while ( iter != m_verticalCurveXPoints.end() )
      {
	if ( m_traits->compare_y(*iter, p) == SMALLER )
	  ++iter; 
	else
	  break;
      }
      if ( iter == m_verticalCurveXPoints.end() )
	m_verticalCurveXPoints.push_back(p);
      else if (!m_traits->point_is_same(p, *iter)) {
	m_verticalCurveXPoints.insert(iter, p);
      }
    }
  }

  /*! 
    Returns a referece to the list of intersection points on the 
    vertical curves passign through the event. If no vertical curves 
    pass through the event or no intersection curves exist, the list 
    will be empty.
    @return a reference to the list of points.
  */
  VerticalXPointList &getVerticalXPointList() {
    return m_verticalCurveXPoints;
  }

  /*! Mark the event as an intersection point at an interior of a curve.
   */
  void markInternalIntersectionPoint() {
    m_isInternalIntersectionPoint = true;
  }

  /*!
    @return returns true if the event is an intersection point at the 
    interior of at least one of the curves passing throuogh the event 
    point.
   */
  bool isInternalIntersectionPoint() const {
    return m_isInternalIntersectionPoint;
  }

  /*! 
    @return true if the any two curves in the event overlap, false otherwise.
  */
  bool doesContainOverlap() const {
    return m_containsOverlap;
  }

#ifndef NDEBUG
  void Print();
  void PrintVerticalXPoints();
#endif
 
private:

  /*! Whenever a new curve is added to the event at the initialization 
    stage, the right most end point to the left of the event point is 
    updated.
    Precondition: the event is either the source or destination of the curve.
    @param curve a pointer to a new curve added to the event.
  */
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

  /*! The point of the event */
  Point_2 m_point;

  /*! A pointer to a traits class */
  Traits *m_traits;

  /*! A list of curves on the left side of the event, sorted by their y value
      to the left of the point */
  SubcurveContainer *m_leftCurves;

  /*! A list of curves on the right side of the event, sorted by their y value
      to the right of the point */
  SubcurveContainer *m_rightCurves;

  /*! The rightmost curve end point that is to the left of the event
      point. This point is used as a reference point when curves are compared
      to the left of the event point. 
  */
  Point_2 m_rightmostPointToLeft;

  /*! An indication whether this event has been initialized. The event is
      initialized after the first curve has been added to the left of the 
      event. 
  */
  bool m_isInitialized;

  /*! a list of vertical curves going through this event */
  VerticalCurveList m_verticalCurves; 

  /*! a list of intersection points on the vertical curves */
  VerticalXPointList m_verticalCurveXPoints;

  /*! a flag that inidcates whether the event is an "interior" intersection 
      point, or just an end point of all curves passing through it.
  */
  bool m_isInternalIntersectionPoint;

  /*! true if any two curves passing through the event overlap. */
  bool m_containsOverlap;



#ifndef NDEBUG
public:
  int id;
#endif
  
};





#ifndef NDEBUG
template<class SweepLineTraits_2, class CurveWrap>
void 
Sweep_line_event<SweepLineTraits_2, CurveWrap>::
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
  for ( SubCurveIter iter1 = m_rightCurves->begin() ;
	iter1 != m_rightCurves->end() ; ++iter1 )
  {
    std::cout << "\t";
    (*iter1)->Print();
    std::cout << "\n";
  }
  std::cout << std::endl;
}

template<class SweepLineTraits_2, class CurveWrap>
void 
Sweep_line_event<SweepLineTraits_2, CurveWrap>::
PrintVerticalXPoints()
{
  std::cout << "Vertical intersection points for " << m_point << ":\n";
  std::list<Point_2>::iterator iter = m_verticalCurveXPoints.begin();
  while ( iter != m_verticalCurveXPoints.end() )
  {
    std::cout << "\t" << *iter << "\n";
    ++iter;
  }
}
 
#endif // NDEBUG

CGAL_END_NAMESPACE

#endif // CGAL_SWEEP_LINE_EVENT_H
