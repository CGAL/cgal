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
// file          : include/CGAL/Sweep_line_tight_2.h
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
#ifndef CGAL_SWEEP_LINE_BASE_H
#define CGAL_SWEEP_LINE_BASE_H

#include <CGAL/Sweep_line_2/Sweep_line_functors.h>

#include <map>
#include <set>


#ifndef VERBOSE

#define SL_DEBUG(a)
#define PRINT_INSERT(a)
#define PRINT_ERASE(a)
#define PRINT_NEW_EVENT(p, e) 
#define DBG(a)
#define STORE_RESULT(a)

#else

#define SL_DEBUG(a) {a}
#define PRINT_INSERT(a) { std::cout << "+++ inserting "; \
                          (a)->Print(); \
                          std::cout << "    currentPos = " << m_currentPos \
                                    << "\n"; \
                          }
#define PRINT_ERASE(a)  { std::cout << "--- erasing " ; \
                          (a)->Print(); }
#define PRINT_NEW_EVENT(p, e) { std::cout << "%%% a new event was created at " \
                                          << (p) << std::endl; \
                                (e)->Print(); }
#define DBG(a) { std::cout << a << std::endl; }
#define STORE_RESULT(a) {a}
#endif


CGAL_BEGIN_NAMESPACE


/*!
  Sweep_line_tight_2 is a class that implements the sweep line algorithm
  based on the algorithm of Bentley and Ottmann.
  It extends the algorithm to support not only segments but polylines and 
  general curves as well.
  The curves are defined by the traits class passed at the template argument.

  The algorithm is also extended to support the following degenerate cases:
  - non x-monotone curves
  - vertical segments
  - moe then two segments intersecting at the same point
  - curves beginning and ending on other curves.

  General flow:
  Init
  For each event
    First pass
    Handle vertical curve (bottom)
    Handle left curves 
    Handle vertical curve (top)
    Handle right curves
  End

  Convensions throuout the code:
  In order to make the code as readable as possible, some convensions were 
  made in regards to variable naming:
    xp - is the intersection point between two curves
    slIter - an iterator to the status line, always points to a curve.

*/



template <class CurveInputIterator,  class SweepLineTraits_2,
          class SweepEvent, class CurveWrap>
class Sweep_line_tight_2
{
public:
  typedef SweepLineTraits_2 Traits;
  typedef typename Traits::Point_2 Point_2;
  typedef typename Traits::Curve_2 Curve_2;
  typedef typename Traits::X_curve_2 X_curve_2;

  typedef SweepEvent Event;
  typedef Point_less_functor<Point_2, Traits> PointLess;
  typedef std::map<const Point_2, Event*, PointLess> EventQueue;
  typedef typename EventQueue::iterator EventQueueIter;
  typedef typename EventQueue::value_type EventQueueValueType;
  typedef std::vector<Event*> EventPtrContainer;
  typedef typename EventPtrContainer::iterator EventPtrContainerIter;

  typedef typename Event::SubCurveIter EventCurveIter;

  typedef CurveWrap Subcurve;
  typedef typename std::list<Subcurve*> SubCurveList;
  typedef typename SubCurveList::iterator SubCurveListIter;

  typedef Status_line_curve_less_functor<Traits, Subcurve> StatusLineCurveLess;
  typedef typename std::set<Subcurve*, StatusLineCurveLess> StatusLine;
  typedef typename StatusLine::iterator StatusLineIter;

  typedef typename Event::VerticalCurveList VerticalCurveList;
  typedef typename Event::VerticalCurveListIter VerticalCurveListIter;
  typedef typename Event::VerticalXPointList VerticalXPointList;
  typedef typename Event::VerticalXPointListIter VerticalXPointListIter;

  typedef std::list<Event *> EventList;
  typedef typename EventList::iterator EventListIter;

  typedef std::list<X_curve_2> CurveList;
  typedef typename CurveList::iterator CurveListIter;

  class  SweepLineGetSubCurves {};
  class  SweepLineGetPoints {};
  class  SweepLineGetInterCurveList {};

  Sweep_line_tight_2() : m_traits(new Traits()), m_traitsOwner(true) {}
  Sweep_line_tight_2(Traits *t) : m_traits(t), m_traitsOwner(false) {}

  ~Sweep_line_tight_2();

  /*!
   *  Given a container of curves, this function returns a list of curves
   *  that are created by intersecting the input curves.
   *  \param curves_begin the input iterator that points to the first curve 
   *                      in the range.
   *  \param curves_end the input past-the-end iterator of the range.
   *  \param subcurves an iterator to the first curve in the range
   *                   of curves created by intersecting the input curves.
   *  \param overlapping indicates whether there are overlapping curves
   *                     in the input range. Defaults to false.
   */
  template <class OutpoutIterator>
  void  get_subcurves(CurveInputIterator begin, CurveInputIterator end, 
		      OutpoutIterator subcurves, bool overlapping = false)
  { 
    SL_DEBUG(std::cout << "*****************************************\n"
                       << "         get_subcurves \n"
	               << "*****************************************\n";)
    Init(begin, end);
    SL_DEBUG(
      PrintSubCurves();
      PrintEventQueue();
    )
    m_overlapping = overlapping;
    Sweep(subcurves, SweepLineGetSubCurves());
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
    SL_DEBUG(std::cout << "*****************************************\n"
                       << "         get_intersection_points \n"
	               << "*****************************************\n";)
    Init(begin, end);
    SL_DEBUG(
      PrintSubCurves();
      PrintEventQueue();
    )
    Sweep(includeEndPoints, points, SweepLineGetPoints());
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
  void  get_intersecting_curves(CurveInputIterator curves_begin, 
				CurveInputIterator curves_end, 
				OutputIterator intersecting_curves,
				bool endpoints = true)
  { 
    SweepLineGetInterCurveList tag;
  }

protected:

  void Init(CurveInputIterator begin, CurveInputIterator end);
  void InitCurve(X_curve_2 &curve);

  /*! The main loop to calculate intersections among the curves
   *  Looping over the events in the queue, for each event we first
   *  handle the curves that are tothe left of the event point (i.e., 
   *  curves that we are done with), and then we look at the curves 
   *  to the right of the point, which means we attept to find intersections
   *  between them and their neighbours on the sweep line.
   */
  template <class OutpoutIterator>
  void Sweep(OutpoutIterator out, SweepLineGetSubCurves tag)
  {
    EventQueueIter eventIter = m_queue->begin();
    m_prevPos = eventIter->first;
    Point_2 referencePoint;

    while ( eventIter != m_queue->end() )
    {
      const Point_2 *p = &(eventIter->first);
      if ( m_traits->compare_x(m_sweepLinePos, *p) == SMALLER )
        m_prevPos = m_sweepLinePos;
      m_sweepLinePos = *p;
      m_currentPos = *p;
      referencePoint = *p;
      m_verticals.clear();
      m_verticalSubCurves.clear();

      while (eventIter != m_queue->end() && 
             m_traits->compare_x(eventIter->first, referencePoint) == EQUAL) {
        p = &(eventIter->first);
        m_currentEvent = eventIter->second;
        SL_DEBUG(std::cout << "------------- " << *p << " --------------"
                           << std::endl;
                 PrintStatusLine();
                 m_currentEvent->Print();
        )

        FirstPass();
	HandleVerticalCurveBottom(tag);
	HandleVerticalOverlapCurves();
        HandleLeftCurves(out, tag);

        m_miniq.push_back(m_currentEvent);
        ++eventIter;
      }
      m_queue->erase(m_queue->begin(), eventIter);

      EventListIter itt = m_miniq.begin();
      while ( itt != m_miniq.end())
      {
        m_currentEvent = *itt;
	HandleVerticalCurveTop(out, tag);
        HandleRightCurves();
        ++itt;
      }
      m_miniq.clear();
      eventIter = m_queue->begin();
    }

    // intersect vertical curves...
    IntersectVerticalCurves();
  }

  void IntersectVerticalCurves()
  {
  }

  /*! The main loop to calculate intersections among the curves
   *  Looping over the events in the queue, for each event we first
   *  handle the curves that are tothe left of the event point (i.e., 
   *  curves that we are done with), and then we look at the curves 
   *  to the right of the point, which means we attept to find intersections
   *  between them and their neighbours on the sweep line.
   */
  template <class OutpoutIterator>
  void Sweep(bool includeEndPoints, OutpoutIterator out,
	     SweepLineGetPoints tag)
  {
    EventQueueIter eventIter = m_queue->begin();
    m_prevPos = eventIter->first;
    Point_2 referencePoint;

    while ( eventIter != m_queue->end() )
    {
      const Point_2 *p = &(eventIter->first);
      if ( m_traits->compare_x(m_sweepLinePos, *p) == SMALLER )
        m_prevPos = m_sweepLinePos;
      m_sweepLinePos = *p;
      m_currentPos = *p;
      referencePoint = *p;
      m_verticals.clear();
      while (eventIter != m_queue->end() && 
             m_traits->compare_x(eventIter->first, referencePoint) == EQUAL) {
        p = &(eventIter->first);
        m_currentEvent = eventIter->second;
        SL_DEBUG(std::cout << "------------- " << *p << " --------------"
                 << std::endl;
                 PrintStatusLine();
                 m_currentEvent->Print();
        )

        FirstPass();
        HandleVerticalCurveBottom(tag);
	HandleVerticalOverlapCurves();
	HandleLeftCurves(includeEndPoints, out, tag);

        m_miniq.push_back(m_currentEvent);
        ++eventIter;
      }

      m_queue->erase(m_queue->begin(), eventIter);

      EventListIter itt = m_miniq.begin();
      while ( itt != m_miniq.end())
      {
        m_currentEvent = *itt;
	HandleVerticalCurveTop(includeEndPoints, out, tag);
        HandleRightCurves();
        ++itt;
      }
      m_miniq.clear();
      eventIter = m_queue->begin();
    }
  }

  void FirstPass();
  void HandleVerticalCurveBottom(SweepLineGetSubCurves &tag);
  void HandleVerticalCurveBottom(SweepLineGetPoints &tag);
  void HandleVerticalOverlapCurves();


  /*!
   *  Handle a vertical curve when the event being processed is the top end 
   *  of the curve. In this situation, the event contains a list of intersection
   *  points on the vertical curve. We go through this list and outpt the 
   *  subcurves induced by these intersection points.
   *  If the curve is not vertical, returns without doing anything.
   * 
   *  @param out an iterator to the output
   *  @param tag a tag that indicates the version of the method
   */
  template <class OutpoutIterator>
  void HandleVerticalCurveTop(OutpoutIterator out, SweepLineGetSubCurves &tag)
  {
    SL_DEBUG(std::cout<<"HandleVerticalCurveTop... (" 
	              << m_currentEvent->getPoint() << ")\n";)
    if ( !m_currentEvent->doesContainVerticalCurve() ) {
      SL_DEBUG(std::cout<<"exiting\n ";)
      return;
    }
    SL_DEBUG(std::cout<<"\n ";)

    VerticalCurveList &vcurves = m_currentEvent->getVerticalCurves();
    VerticalCurveListIter vciter = vcurves.begin();

    while ( vciter !=vcurves.end() )
    {

      Subcurve *vcurve = *vciter;
      const Point_2 &topPoint = m_currentEvent->getPoint();
      // if this is the bottom point, nothing to do here
      if ( vcurve->isBottomEnd(topPoint)) {
	SL_DEBUG(std::cout<<"this is the bottom. skipping.\n";)
	++vciter;
	continue;
      }

      SL_DEBUG(std::cout<<"handling top point of vertical curve\n";)


      // the following while loop comes to handle  | 
      // in the case where a new curve begins at   |------
      // a vertical curve                          |

      // find the "position" of the curve of the status line
      StatusLineIter slIter = m_statusLine->lower_bound(vcurve);
      
      if ( slIter != m_statusLine->end() ) 
      {
	SL_DEBUG(std::cout<<"starting at curve \n";)
	SL_DEBUG((*slIter)->Print();)

	while ( slIter != m_statusLine->end() &&
		m_traits->curve_get_point_status((*slIter)->getCurve(), 
						 topPoint) 
		== Traits::ABOVE_CURVE &&
		m_traits->curve_get_point_status((*slIter)->getCurve(), 
						 vcurve->getBottomEnd()) 
		== Traits::UNDER_CURVE )
	{
	  SL_DEBUG(std::cout<<"checking \n";)
	  SL_DEBUG((*slIter)->Print();) 
	  if ( m_traits->compare_x((*slIter)->getLeftEnd(), topPoint) == EQUAL)
	  {
	    m_currentEvent->addVerticalCurveXPoint((*slIter)->getLeftEnd(), true);
	  }
	  ++slIter;
	}   
      }

      // now we go over the list of intersection points on the vertical
      // curve in at the event and process them...
      SL_DEBUG(std::cout<<"handling the splitting now\n";)
      VerticalXPointList &pointList = m_currentEvent->getVerticalXPointList();
      if ( pointList.empty() )
      {
	AddVerticalCurveToOutput(out, vcurve->getCurve());
	//*out = vcurve->getCurve();
	//++out;
	++vciter;
	continue;
      }
    
      X_curve_2 a, b, c;
      a = vcurve->getCurve();
      SL_DEBUG(std::cout << "there are " << pointList.size() << " points\n";)
      SL_DEBUG(m_currentEvent->PrintVerticalXPoints();)
      for ( VerticalXPointListIter i = pointList.begin() ;
	    i != pointList.end(); ++i )
      {
	SL_DEBUG(std::cout<< "splitting: " << a << " at " << *i ;)
	if ( !vcurve->isPointInRange(*i) )
	{
	  SL_DEBUG(std::cout << " not !\n";)
	  continue;
	}
	SL_DEBUG(std::cout << " yes! \n";)
	m_traits->curve_split(a, b, c, *i);
	if ( vcurve->isSourceLeftToTarget()) {
	  AddVerticalCurveToOutput(out, b);
	  //*out = b; ++out;
	  a = c;
	} else {
	  AddVerticalCurveToOutput(out, c);
	  //*out = c; ++out;
	  a = b;
	}
      }
      if ( vcurve->isSourceLeftToTarget() ) {
	AddVerticalCurveToOutput(out, c);
	//*out = c; ++out;
      }
      else {
	AddVerticalCurveToOutput(out, b);
	//*out = b; ++out;
      }
      ++vciter;
    }
  }

  /*! Adds a new curve to the output list. If the overlapping flag is false,
      each unique curve is reported only once.

      @param out an iterator to the output container
      @param cv the curve to be added
  */
  template <class OutpoutIterator>
  void AddCurveToOutput(const X_curve_2 &cv, Subcurve *curve, 
			OutpoutIterator out)
  {
    static Subcurve *prevCurve = 0;
    static X_curve_2 prevXCv;

    if ( m_overlapping ) {
      *out = cv;
      ++out;
    } else {
      if ( prevCurve && SimilarCurves(cv, prevXCv)) {
	SL_DEBUG(std::cout << " curve already reported... " << std::endl;)
	return;
      }
      prevCurve = curve;
      prevXCv = cv;
      *out = cv;
      ++out;
    }
  }

  template <class OutpoutIterator>
  void AddVerticalCurveToOutput(OutpoutIterator out, 
				const X_curve_2 &cv)
  {
    if ( m_overlapping ) {
      *out = cv;
      ++out;
    } else {
      if ( VerticalSubCurveExists(cv)) {
	SL_DEBUG(std::cout << " curve already reported... " << std::endl;)
	return;
      }
      m_verticalSubCurves.push_back(cv);
      *out = cv;
      ++out;
    }
  }

  /*! 
   * Returns true if the point is in the interior of the curve.
   * Returns false if the point is outside the range of the curve or
   * if the point is either the source or the target of the curve.
   * @return true if the point is int he interior of the curve.
   */
  bool isPointInCurveInterior(const X_curve_2 &c, const Point_2 &p)
  {
    if ( m_traits->curve_get_point_status(c, p) != Traits::ON_CURVE )
      return false;
    if ( isEndPoint(p) )
      return false;
    return true;
  }

  /*!
   *  Handle a vertical curve when the event being processed is the top end 
   *  of the curve. In this situation, the event contains a list of intersection
   *  points on the vertical curve. We go through this list and output the 
   *  intersection points.
   *  If the curve is not vertical, returns without doing anything.
   * 
   *  @param out an iterator to the output
   *  @param tag a tag that indicates the version of the method
   */
  template <class OutpoutIterator>
  void HandleVerticalCurveTop(bool includeEndPoints, OutpoutIterator out, 
			      SweepLineGetPoints &tag)
  {
    SL_DEBUG(std::cout<<"HandleVerticalCurveTop... ";)
    if ( !m_currentEvent->doesContainVerticalCurve() )
    {
      SL_DEBUG(std::cout<<"exiting\n ";)
      return;
    }
    SL_DEBUG(std::cout<<"\n ";)


    VerticalCurveList &vcurves = m_currentEvent->getVerticalCurves();
    VerticalCurveListIter vciter = vcurves.begin();

    while ( vciter != vcurves.end() )
    {

      Subcurve *vcurve = *vciter; 
      const Point_2 &topPoint = m_currentEvent->getPoint();
      if ( vcurve->isBottomEnd(topPoint)) {
	SL_DEBUG(std::cout<<"this is the bottom. skipping.\n";)
	++vciter;
	continue;
      }

      SL_DEBUG(std::cout<<"handling top point of vertical curve\n";)
      StatusLineIter slIter = m_statusLine->lower_bound(vcurve);
      if ( slIter != m_statusLine->end() )
      {
	SL_DEBUG(std::cout<<"starting at curve \n";)
	SL_DEBUG((*slIter)->Print();)
	  
        const Point_2 &bottomPoint = vcurve->getBottomEnd();

	while ( slIter != m_statusLine->end() &&
		m_traits->curve_get_point_status((*slIter)->getCurve(), topPoint) 
		== Traits::ABOVE_CURVE &&
		m_traits->curve_get_point_status((*slIter)->getCurve(), 
						 bottomPoint) 
		== Traits::UNDER_CURVE )
	{
	  SL_DEBUG(std::cout<<"checking \n";)
	    SL_DEBUG((*slIter)->Print();) 
	    if ( m_traits->compare_x((*slIter)->getLeftEnd(),topPoint) == EQUAL)
	    {
	      m_currentEvent->addVerticalCurveXPoint((*slIter)->getLeftEnd());
	      // if this point was not an event, we need to report the point
	      // test40/42
	      if ( !includeEndPoints && 
		   !isInternalXPoint((*slIter)->getLeftEnd()))
	      {
		*out = (*slIter)->getLeftEnd(); ++out;
	      }
	    }
	  ++slIter;
	}   
      }
      ++vciter;
    }
  }

  /*! For each left-curve, if it is the "last" subcurve, i.e., the 
   * event point is the right-edge of the original curve, the 
   * last sub curve is created and added to the result. Otherwise
   * the curve is added as is to the result.
   */
  template <class OutpoutIterator>
  void HandleLeftCurves(OutpoutIterator out, 
			SweepLineGetSubCurves &tag)
  {
    SL_DEBUG(std::cout << "Handling left curve" << std::endl;)
    SL_DEBUG(m_currentEvent->Print();)
    EventCurveIter leftCurveIter = m_currentEvent->leftCurvesBegin();
    m_currentPos = m_prevPos;
    const Point_2 &eventPoint = m_currentEvent->getPoint();

    while ( leftCurveIter != m_currentEvent->leftCurvesEnd() )  // ** fix here
    {
      Subcurve *leftCurve = *leftCurveIter; 
      const X_curve_2 &cv = leftCurve->getCurve();
      const Point_2 &lastPoint = leftCurve->getLastPoint();

      if ( leftCurve->isSource(eventPoint))
      {
        if ( !leftCurve->isTarget(lastPoint) )
        {
          X_curve_2 a,b;
          m_traits->curve_split(cv, a, b, lastPoint);
          AddCurveToOutput(a, leftCurve, out);
        } else {
          AddCurveToOutput(cv, leftCurve, out);
        }
      } else if ( leftCurve->isTarget(eventPoint))
      {
        if ( !leftCurve->isSource(lastPoint))
        {
          X_curve_2 a,b;
          m_traits->curve_split(cv, a, b, lastPoint);
          AddCurveToOutput(b, leftCurve, out);
        } else {
          AddCurveToOutput(cv, leftCurve, out);
        }

      } else { 
        X_curve_2 a,b;
        if ( leftCurve->isSource(lastPoint)) {
          m_traits->curve_split(cv, a, b, eventPoint);
          AddCurveToOutput(a, leftCurve, out);
        } else if ( leftCurve->isTarget(lastPoint)) {
          m_traits->curve_split(cv, b, a, eventPoint);
          AddCurveToOutput(a, leftCurve, out);
        } else {
          const X_curve_2 &lastCurve = leftCurve->getLastCurve();
          if ( leftCurve->isSourceLeftToTarget() ) {
            m_traits->curve_split(lastCurve, a, b, eventPoint);
            AddCurveToOutput(a, leftCurve, out);
          } else {
            m_traits->curve_split(lastCurve, b, a, eventPoint);
            AddCurveToOutput(a, leftCurve, out);
          }
        }
        leftCurve->setLastPoint(eventPoint);
        leftCurve->setLastCurve(b); 
      }

      // before deleting check new neighbors that will become after deletion
      StatusLineIter sliter = 
	IntersectNeighboursAfterRemoval(leftCurve);

      m_currentPos = m_prevPos;
      m_statusLine->erase(sliter);
      ++leftCurveIter;
    }
  }


  /*! For each left-curve, if it is the "last" subcurve, i.e., the 
   * event point is the right-edge of the original curve, the 
   * last sub curve is created and added to the result. Otherwise
   * the curve is added as is to the result.
   */
  template <class OutpoutIterator>
  void HandleLeftCurves(bool includeEndPoints,
			OutpoutIterator out, SweepLineGetPoints &tag)
  {
    SL_DEBUG(std::cout << "Handling left curve" << std::endl;)
    SL_DEBUG(m_currentEvent->Print();)
    const Point_2 &eventPoint = m_currentEvent->getPoint();
    if ( !m_currentEvent->hasLeftCurves() )
    {
      if ( includeEndPoints || m_currentEvent->isInternalIntersectionPoint()) {
        *out = eventPoint; ++out;    
      }
      return;
    }

    // delete the curve from the status line
    EventCurveIter leftCurveIter = m_currentEvent->leftCurvesBegin();
    m_currentPos = m_prevPos;
    while ( leftCurveIter != m_currentEvent->leftCurvesEnd() )
    {
      // before deleting check new neighbors that will become after deletion
      StatusLineIter sliter = 
	IntersectNeighboursAfterRemoval(*leftCurveIter);

      PRINT_ERASE((*leftCurveIter));
      m_currentPos = m_prevPos;
      m_statusLine->erase(sliter);
      ++leftCurveIter;
    }

    if ( includeEndPoints || m_currentEvent->isInternalIntersectionPoint() )
    {
      *out = eventPoint; ++out;
    }
  }

  

  void HandleRightCurves();
  bool  Intersect(Subcurve *c1, Subcurve *c2);
  void IntersectCurveGroup(Subcurve *c1, SubCurveList &mylist);

  bool isInternalXPoint(const Point_2 &p);
  bool HandleVerticalCurveXAtEnd(Subcurve *vcurve, Subcurve *curve, 
				 Event *topEndEvent, SweepLineGetSubCurves tag);
  bool HandleVerticalCurveXAtEnd(Subcurve *vcurve, Subcurve *curve, 
				 Event *topEndEvent, SweepLineGetPoints tag);

  /*!
   * When a curve is removed from the status line for good, its top and
   * bottom neighbors become neighbors. This method finds these cases and
   * looks for the itnersection point, if one exists.
   * @param leftCurve a pointer to the curve that is about to be deleted
   * @return an iterator to the position where the curve will be removed from.
   */

  StatusLineIter IntersectNeighboursAfterRemoval(Subcurve *leftCurve)
  {
    SL_DEBUG(PrintStatusLine();)
    SL_DEBUG(leftCurve->Print();)

    StatusLineIter sliter = m_statusLine->find(leftCurve);
    if ( !leftCurve->isEndPoint(m_currentEvent->getPoint()))
      return sliter;
  
    m_currentPos = m_prevPos;
    assert(sliter!=m_statusLine->end());
    StatusLineIter end = m_statusLine->end(); --end;
    if ( sliter != m_statusLine->begin() && sliter != end ) 
    {
      SubCurveList mylist;
      StatusLineIter prev = sliter; --prev;
    
      // collect all curves that overlap with *prev
      StatusLineIter tmp = prev;
      mylist.push_back(*prev);
      while ( tmp != m_statusLine->begin() ) 
      {
	--tmp;
	if ( DoCurvesOverlap(*prev, *tmp))
	  mylist.push_back(*tmp);
	else
	  break;
      }
    
      StatusLineIter next = sliter; ++next;
    
      // intersect *next with the the *prev curve and all overlaps
      tmp = next;
      IntersectCurveGroup(*tmp, mylist);
      
      // if there are curves that overlap with the *next curve, intersect
      // them with the *prev curve and all overlaps
      ++tmp;
      while ( tmp != m_statusLine->end() ) 
      {
	if ( DoCurvesOverlap(*next, *tmp))
	{
	  IntersectCurveGroup(*tmp, mylist);
	  ++tmp;
	}
	else
	  break;
      }
    }
    
    return sliter;
  } 

  bool DoCurvesOverlap(Subcurve *c1, Subcurve *c2);
  bool SimilarCurves(const X_curve_2 &a, const X_curve_2 &b);
  bool VerticalSubCurveExists(const X_curve_2 &a);

  void PrintEventQueue();
  void PrintSubCurves();
  void PrintStatusLine();
  void PrintVerticals();

protected:
  /*! a pointer to a traits object */
  Traits *m_traits;

  /*! an indication to whether the traits should be deleted in the distructor */
  bool m_traitsOwner;

  /*! if false, overlapping subcurves are reported only one. 
    Otherwise, they are reported as many times as they appeard. */
  bool m_overlapping;

  /*! used to hold all event, processed or not, so that they can 
`     all be deleted at the end */
  EventPtrContainer m_events; 

  /*! the queue of events (intersection points) to handle */
  EventQueue *m_queue;

  /*! The subcurves, as created on the fly */
  SubCurveList m_subCurves;

  /*! The status line */
  StatusLine *m_statusLine;

  /*! A reference point that is used for comapring curves. It is used
      when inserting/erasing curves from the status line. */
  Point_2 m_currentPos;

  /*! A reference point that is used for comapring curves */
  Point_2 m_prevPos;

  /*! The current position (in X) of the sweep line */
  Point_2 m_sweepLinePos;

  /*! if non x-monotone are specified, this hold the x-monotone 
    curves created when splitting them into x-monotone curves. */
  std::vector<X_curve_2> m_xcurves;

  /*! a pointer to thecurrent event */
  Event *m_currentEvent;

  /*! a queue that holds all the events that have the same x coordinate as 
      the status line. */
  EventList m_miniq;

  /*! a list of vertical curves at the x coordinate of the current event 
      point.*/
  SubCurveList m_verticals;
  CurveList m_verticalSubCurves;

  /*! a counter the is used to assign unique ids to the curves. */
  int m_curveId;

#ifndef NDEBUG
  int m_eventId;
#endif
};

template <class CurveInputIterator,  class SweepLineTraits_2,
         class SweepEvent, class CurveWrap>
Sweep_line_tight_2<CurveInputIterator,SweepLineTraits_2,SweepEvent,CurveWrap>::
~Sweep_line_tight_2() 
{
  if ( m_traitsOwner ) delete m_traits;

  for ( SubCurveListIter sci = m_subCurves.begin() ; 
	sci != m_subCurves.end() ; ++sci)
  {
    delete *sci;
  }

  for ( EventPtrContainerIter ei = m_events.begin();
	ei != m_events.end() ; ++ei)
  {
    delete *ei;
  }
  delete m_queue;
  delete m_statusLine;
}



/*! Initializes the data structures to work with:
 *  - x-monotonize the inf\put curves
 *  - for each end point of each curve create an event
 *  - initialize the event queue
 *  -
 */
template <class CurveInputIterator,  class SweepLineTraits_2,
         class SweepEvent, class CurveWrap>
inline void 
Sweep_line_tight_2<CurveInputIterator, SweepLineTraits_2,SweepEvent,CurveWrap>::
Init(CurveInputIterator begin, CurveInputIterator end)
{
  PointLess pred(m_traits);
  m_queue = new EventQueue(pred);
  StatusLineCurveLess slcurveless(m_traits);
  m_statusLine = new StatusLine(slcurveless);
#ifndef NDEBUG
  m_eventId = 0;
#endif
  m_curveId = 0;

  int count = 0;
  CurveInputIterator iter;
  for ( iter = begin ; iter != end ; ++iter)
  {
    if ( m_traits->is_x_monotone(*iter) ) 
      InitCurve(*iter);
    else
    {
      std::list<X_curve_2> xcurves;
      m_traits->make_x_monotone(*iter, xcurves);
      SL_DEBUG(
      std::cout << "curve " << *iter << " was split into " 
                << xcurves.size() << " curves." << std::endl;
      )

      for ( typename std::list<X_curve_2>::iterator i = xcurves.begin();
	    i != xcurves.end() ; ++i )
      {
	m_xcurves.push_back(*i);
	InitCurve(m_xcurves[count]);
	count++;
      }
    }
  }
}



/*! Given an x-monotone curve, create events for each end (if 
 *  one doesn't exist already). 
 *  For each curve create a Subcurve instance.
 */
template <class CurveInputIterator,  class SweepLineTraits_2,
         class SweepEvent, class CurveWrap>
inline void 
Sweep_line_tight_2<CurveInputIterator, SweepLineTraits_2,SweepEvent,CurveWrap>::
InitCurve(X_curve_2 &curve)
{
  const Point_2 &source = m_traits->curve_source(curve);
  const Point_2 &target = m_traits->curve_target(curve);
  Event *e = 0;
  
  Subcurve *subCv = new Subcurve(m_curveId++, curve, &m_currentPos, m_traits);
  m_subCurves.push_back(subCv);
  
  // handle the source point
  EventQueueIter eventIter = m_queue->find(source);
  if ( eventIter != m_queue->end() ) {
    SL_DEBUG(std::cout << "event " << source << " already exists\n";)
    e = eventIter->second;
  } else  {
    e = new Event(source, m_traits); 
#ifndef NDEBUG
    e->id = m_eventId++;
#endif
    m_events.push_back(e);
    m_queue->insert(EventQueueValueType(source, e));
  }
  e->addCurve(subCv);
  PRINT_NEW_EVENT(source, e);
    
  // handle the target point
  eventIter = m_queue->find(target);
  if ( eventIter != m_queue->end() ) {
    SL_DEBUG(std::cout << "event " << target << " already exists\n";)
    e = eventIter->second;
  } else  {
    e = new Event(target, m_traits); 
#ifndef NDEBUG
    e->id = m_eventId++;
#endif
    m_events.push_back(e);
    m_queue->insert(EventQueueValueType(target, e));
  }
  e->addCurve(subCv);
  PRINT_NEW_EVENT(target, e);
}



/*! This pass comes to take care of cases in which we reach an
 *  event point that is the left end of a curve that starts on another curve.
 *  For example:
 *       /      ------         /
 *      /----      \          /
 *     /            \     ----------
 *  This method is called before the left curves of an event are handled.
 */
template <class CurveInputIterator,  class SweepLineTraits_2,
         class SweepEvent, class CurveWrap>
inline void 
Sweep_line_tight_2<CurveInputIterator, SweepLineTraits_2,SweepEvent,CurveWrap>::
FirstPass()
{
  if ( m_statusLine->size() == 0 )
    return;

  const Point_2& p = m_currentEvent->getPoint();

  DBG("First pass");
  EventCurveIter rightIter = m_currentEvent->rightCurvesBegin();
  m_currentPos = m_sweepLinePos;
  while ( rightIter != m_currentEvent->rightCurvesEnd())
  {
    // if the point is not an end point, skip it
    if ( !(*rightIter)->isEndPoint(p)) {
      ++rightIter;
      continue;
    }

    StatusLineIter slIter = m_statusLine->lower_bound(*rightIter);
    StatusLineIter prev = slIter;
    StatusLineIter next = slIter;
 
    // check the curve that comes before the rightCurve on the status line.
    // check also all of the overlaps...
    if ( slIter != m_statusLine->begin() )
    { 
      --prev;
      while ( m_traits->curve_get_point_status((*prev)->getCurve(), p) ==
	      Traits::ON_CURVE && !(*prev)->isEndPoint(p))
      {
	m_currentEvent->addCurveToRight(*prev);
	m_currentEvent->addCurveToLeft(*prev, m_prevPos);
	if ( prev == m_statusLine->begin() )
	  break;
	--prev;
      }
    }

    // check the curve that comes after the rightCurve on the status line.
    // check also all of the overlaps...
   if ( slIter != m_statusLine->end() )
    {
      while ( m_traits->curve_get_point_status((*next)->getCurve(), p) ==
	      Traits::ON_CURVE && !(*next)->isEndPoint(p))
      {
	m_currentEvent->addCurveToRight(*next);
	m_currentEvent->addCurveToLeft(*next, m_prevPos);
	++next;
	if ( next ==  m_statusLine->end() )
	  break;
      }    
    } 
    ++rightIter;
  }
  SL_DEBUG(m_currentEvent->Print();)
  SL_DEBUG(std::cout << "First pass - done\n" ;)
}



/*!
 * Handles the degenerate case of vertical curves. Most of the cases
 * that occur with vertical curves are handled by this method and 
 * HandleVerticalCurveTop method.
 * When the current event is the bottom end of a vertical curve, we look
 * for intersection points between the vertical curve and any curve
 * in the status line that in the y-range that is defined by the bottom 
 * and top ends of the vertical curve. When those are found, we create
 * new events, unless ones already exist, in which case we update the events.
 * 
 * @param tag a tag that indicates the version of this method
 * \sa HandleVerticalCurveTop
 */
template <class CurveInputIterator,  class SweepLineTraits_2,
         class SweepEvent, class CurveWrap>
inline void 
Sweep_line_tight_2<CurveInputIterator, SweepLineTraits_2,SweepEvent,CurveWrap>::
HandleVerticalCurveBottom(SweepLineGetSubCurves &tag)
{
  SL_DEBUG(std::cout<<"\nHandleVerticalCurveBottom... ("
                    << m_currentEvent->getPoint() << ")\n";)
  if ( !m_currentEvent->doesContainVerticalCurve() )
  {
    SL_DEBUG(std::cout<<" - not vertical - exiting\n ";)
    return;
  }
  SL_DEBUG(std::cout<<"\n ";)

  VerticalCurveList &vcurves = m_currentEvent->getVerticalCurves();
  VerticalCurveListIter vciter = vcurves.begin();
  const Point_2 &currentPoint = m_currentEvent->getPoint();

  SL_DEBUG(std::cout << vcurves.size() << " vertical curves in event\n";)
  while ( vciter != vcurves.end() )
  {
    Subcurve *vcurve = *vciter;
    SL_DEBUG(std::cout << "working on " << vcurve->getCurve() << "\n";)
    if ( vcurve->isTopEnd(currentPoint))
    {
      vciter++;
      continue;
    }
    
    SL_DEBUG(std::cout<<"handling bottom point of vertical curve\n";)
    StatusLineIter slIter = m_statusLine->lower_bound(vcurve);
    if ( slIter == m_statusLine->end() ) {
      SL_DEBUG(std::cout<<"no curves intersecting. exiting\n";)
      vciter++;
      continue;
    }    
    
    SL_DEBUG(std::cout<<"starting at curve \n";)
    SL_DEBUG((*slIter)->Print();)
    const Point_2 &topEnd = vcurve->getTopEnd();
    EventQueueIter topEndEventIter = m_queue->find(topEnd);
    assert(topEndEventIter!=m_queue->end());
    Event *topEndEvent = topEndEventIter->second;

    bool lastEventCreatedHere = false;
    Event *prevEvent = 0;

    while ( slIter != m_statusLine->end() &&
	    m_traits->curve_get_point_status((*slIter)->getCurve(), topEnd) 
	    != Traits::UNDER_CURVE &&
	    m_traits->curve_get_point_status((*slIter)->getCurve(), 
					     currentPoint) 
	    != Traits::ABOVE_CURVE )
    {
      SL_DEBUG(std::cout<<"intersecting with \n";)
      SL_DEBUG((*slIter)->Print();) 
	
      if ( HandleVerticalCurveXAtEnd(vcurve, *slIter, topEndEvent, tag))
      {
	++slIter;
	continue;
      }
      
      // handle a curve that goes through the interior of the vertical curve
      const X_curve_2 &cv1 = vcurve->getCurve();
      const X_curve_2 &cv2 = (*slIter)->getCurve();
      Point_2 p;
      bool res =
	m_traits->nearest_intersection_to_right(cv1, cv2, currentPoint, p, p);
      SL_DEBUG(assert(res==true);)
      res = 0;
      
      EventQueueIter eqi = m_queue->find(p);
      Event *e = 0;
      if ( eqi == m_queue->end() )
      {
	e = new Event(p, m_traits); 
#ifndef NDEBUG
	e->id = m_eventId++;
#endif
	m_events.push_back(e);
	
	e->addCurveToLeft(*slIter, m_sweepLinePos);
	e->addCurveToRight(*slIter);
	PRINT_NEW_EVENT(p, e);
	m_queue->insert(EventQueueValueType(p, e));

	lastEventCreatedHere = true;

      } else {
	e = eqi->second;
	
	// the only time we need to update the event is when the event
	// is created here (which also includes overlapping curves)
	if ( e == prevEvent ) {
	  if ( lastEventCreatedHere )
	  {
	    if ( !(*slIter)->isLeftEnd(p) ) 
	      e->addCurveToLeft(*slIter, m_sweepLinePos);
	    if ( !(*slIter)->isRightEnd(p) ) 
	      e->addCurveToRight(*slIter);
	  } 
	}
	else
	  lastEventCreatedHere = false;

	SL_DEBUG(std::cout << "Updating event \n";)
	SL_DEBUG(e->Print();)
      }
      
      topEndEvent->addVerticalCurveXPoint(p);
      ++slIter;
      prevEvent = e;
    }    
    vciter++;
  }

  SL_DEBUG(std::cout<<"Done Handling vertical\n";)
}



/*!
 * Handles the degenerate case of vertical curves. Most of the cases
 * that occur with vertical curves are handled by this method and the
 * HandleVerticalCurveTop method.
 *
 * When the current event is the bottom end of a vertical curve, we look
 * for intersection points between the vertical curve and any curve
 * in the status line that in the y-range that is defined by the bottom 
 * and top ends of the vertical curve. When those are found, we create
 * new events, unless ones already exist, in which case we update the events.
 * 
 * @param tag a tag that indicates the version of this method
 * \sa HandleVerticalCurveTop
 */
template <class CurveInputIterator,  class SweepLineTraits_2,
         class SweepEvent, class CurveWrap>
inline void 
Sweep_line_tight_2<CurveInputIterator, SweepLineTraits_2,SweepEvent,CurveWrap>::
HandleVerticalCurveBottom(SweepLineGetPoints &tag)
{
  SL_DEBUG(std::cout<<"HandleVerticalCurveBottom... ";)
  if ( !m_currentEvent->doesContainVerticalCurve() )
  {
    SL_DEBUG(std::cout<<"exiting\n ";)
    return;
  }
  SL_DEBUG(std::cout<<"\n ";)

  VerticalCurveList &vcurves = m_currentEvent->getVerticalCurves();
  VerticalCurveListIter vciter = vcurves.begin();
  const Point_2 &currentPoint = m_currentEvent->getPoint();

  while ( vciter != vcurves.end() )
  {
    Subcurve *vcurve = *vciter;
    if ( vcurve->isTopEnd(currentPoint))
    {
      ++vciter;
      continue;
    }

    SL_DEBUG(std::cout<<"handling bottom point of vertical curve\n";);
    StatusLineIter slIter = m_statusLine->lower_bound(vcurve);
    if ( slIter == m_statusLine->end() ) {
      SL_DEBUG(std::cout<<"no curves intersecting. exiting\n";);
      ++vciter;
      continue;
    }    

    SL_DEBUG(std::cout<<"starting at curve \n";);
    SL_DEBUG((*slIter)->Print(););

    const Point_2 &topEnd = vcurve->getTopEnd();
    EventQueueIter topEndEventIter = m_queue->find(topEnd);
    assert(topEndEventIter!=m_queue->end());
    Event *topEndEvent = topEndEventIter->second;

    while ( slIter != m_statusLine->end() &&
	    m_traits->curve_get_point_status((*slIter)->getCurve(), topEnd) 
	    != Traits::UNDER_CURVE &&
	    m_traits->curve_get_point_status((*slIter)->getCurve(), 
					     currentPoint) 
	    != Traits::ABOVE_CURVE )
    {
      SL_DEBUG(std::cout<<"intersecting with \n";)
      SL_DEBUG((*slIter)->Print();) 
	
      if ( HandleVerticalCurveXAtEnd(vcurve, *slIter, topEndEvent, tag))
      {
	++slIter;
	continue;
      }

      Point_2 xp;
      bool res = 
	m_traits->nearest_intersection_to_right(vcurve->getCurve(), 
						(*slIter)->getCurve(), 
						currentPoint, 
						xp, xp);
      SL_DEBUG(assert(res==true);)
      res = 0;
      EventQueueIter eqi = m_queue->find(xp);
      Event *e = 0;
      if ( eqi == m_queue->end() )
      {
	e = new Event(xp, m_traits); 
#ifndef NDEBUG
	e->id = m_eventId++;
#endif
	m_events.push_back(e);
      
	e->addCurveToLeft(*slIter, m_sweepLinePos);
	e->addCurveToRight(*slIter);
	
	PRINT_NEW_EVENT(xp, e);
	m_queue->insert(EventQueueValueType(xp, e));
      } else {
	e = eqi->second;
	e->markInternalIntersectionPoint();
	SL_DEBUG(std::cout << "Updating event \n";)
	SL_DEBUG(e->Print();)
	e->addCurve(vcurve); // test41
      }
      
      topEndEvent->addVerticalCurveXPoint(xp);
      ++slIter;
    }    
    ++vciter;
  }

  SL_DEBUG(std::cout<<"Done Handling vertical\n";)
}



/*!
 * Handles overlapping vertical curves. 
 * If the current event point does not contain vertical curves, nothing is done 
 * here.
 * Fo the current event point, we go through the list of vertical curves 
 * defined in the same x coordinate (m_verticals). For each curve, we check 
 * if the event point is in the interior of the vertical curve. If so, 
 * the event is set to be an intersection point (between the two 
 * vertical curves). 
 * While going through the vertical curves, if we reach a curve that the 
 * event point is above the curve, we remove the curve from the list.
 * 
 * Finally, we go thorugh the vertical curves of the event. If the event 
 * point is the bottom end of a vertical curve, we add the vertical curve 
 * to the list of vertical curves (m_verticals).
 */
template <class CurveInputIterator,  class SweepLineTraits_2,
         class SweepEvent, class CurveWrap>
inline void 
Sweep_line_tight_2<CurveInputIterator, SweepLineTraits_2,SweepEvent,CurveWrap>::
HandleVerticalOverlapCurves()
{
  SL_DEBUG(std::cout<<"\nHandleVerticalOverlapCurves... (" 
                    << m_currentEvent->getPoint() << ")";)

  if ( !m_currentEvent->doesContainVerticalCurve() ) {
    SL_DEBUG(std::cout << "no vertical - exiting\n";)
    return;
  }
  SL_DEBUG(std::cout << "\n";)
  SL_DEBUG(PrintVerticals();)

  const Point_2 &point = m_currentEvent->getPoint();
  SubCurveListIter iter = m_verticals.begin();
  while ( iter != m_verticals.end() )
  {
    Subcurve *curve = *iter;
    typename Traits::Curve_point_status pstatus = 
      m_traits->curve_get_point_status(curve->getCurve(), point);

    if ( pstatus == Traits::ABOVE_CURVE ) {
      iter = m_verticals.erase(iter);

    } else if (!curve->isEndPoint(point)) {
      EventQueueIter eventIter = m_queue->find(curve->getTopEnd());
      assert(eventIter!=m_queue->end());
      (eventIter->second)->addVerticalCurveXPoint(point, true);
      m_currentEvent->markInternalIntersectionPoint();
      ++iter;
    } else {
      ++iter;
    }
  }

  VerticalCurveList &vcurves = m_currentEvent->getVerticalCurves();
  VerticalCurveListIter vciter = vcurves.begin();
  while ( vciter != vcurves.end() )
  {
    Subcurve *vcurve = *vciter;
    if ( vcurve->isBottomEnd(point) ) {
      m_verticals.push_back(vcurve);
    }
    ++vciter;
  }
}



/*! Loop over the curves to the right of the status line and handle them:
 * - if we are at the beginning of the curve, we insert it to the status 
 *   line, then we look if it intersects any of its neighbours.
 * - if we are at an intersection point between two curves, we add them
 *   to the status line and attempt to intersect them with their neighbours. 
 * - We also check to see if the two intersect again to the right of the point.
 */
template <class CurveInputIterator,  class SweepLineTraits_2,
         class SweepEvent, class CurveWrap>
inline void 
Sweep_line_tight_2<CurveInputIterator, SweepLineTraits_2,SweepEvent,CurveWrap>::
HandleRightCurves()
{
  SL_DEBUG(std::cout << "Handling right curves (" ;)
  SL_DEBUG(std::cout << m_currentEvent->getPoint() << ")\n";)
  int numRightCurves = m_currentEvent->getNumRightCurves();
  if ( numRightCurves == 0 )
    return;

  m_currentPos = m_sweepLinePos;
  if ( numRightCurves == 1 )
  {
    SL_DEBUG(std::cout << " - beginning of curve " << std::endl;)

    SL_DEBUG(
      Subcurve *tmp1 = *(m_currentEvent->rightCurvesBegin());
      PRINT_INSERT(tmp1);
    )

    std::pair<StatusLineIter, bool> tmp =
      m_statusLine->insert(*(m_currentEvent->rightCurvesBegin()));
    StatusLineIter slIter = tmp.first;

    SL_DEBUG(PrintStatusLine();)
    // if this is the only curve on the status line, nothing else to do
    if ( m_statusLine->size() == 1 )
      return;

    StatusLineIter prev = slIter;
    StatusLineIter next = slIter;
    ++next;

    SubCurveList mylist;
    if ( slIter != m_statusLine->begin() )
    {
      --prev;
      StatusLineIter tmp = prev;
      mylist.push_back(*prev);
      while ( tmp != m_statusLine->begin() ) 
      {
	--tmp;
	if ( DoCurvesOverlap(*prev, *tmp) )
	  mylist.push_back(*tmp);
	else
	  break;
      }
    }

    if ( next != m_statusLine->end() )
    { 
      StatusLineIter tmp = next;
      mylist.push_back(*next);
      ++tmp;
      while ( tmp != m_statusLine->end() ) 
      {
	if ( DoCurvesOverlap(*next, *tmp) )
	{
	  mylist.push_back(*tmp);
	  ++tmp;
	}
	else
	  break;
      }
    }
    IntersectCurveGroup(*(m_currentEvent->rightCurvesBegin()), mylist);

 
  } else
  {
    SubCurveList mylist;
    SubCurveList prevlist;
    SubCurveList currentlist;


    SL_DEBUG(std::cout << " - intersection point " << std::endl;)
    EventCurveIter firstOne = m_currentEvent->rightCurvesBegin();
    EventCurveIter lastOne = m_currentEvent->rightCurvesEnd(); --lastOne;
    EventCurveIter rightCurveEnd = m_currentEvent->rightCurvesEnd();

    PRINT_INSERT(*firstOne);
    std::pair<StatusLineIter, bool> tmp = m_statusLine->insert(*firstOne);
    StatusLineIter slIter = tmp.first;

    SL_DEBUG(PrintStatusLine();)
    if ( slIter != m_statusLine->begin() )
    { 
      StatusLineIter prev = slIter; --prev;

      // find all curves that are overlapping with the prev curve
      StatusLineIter tmp = prev;
      prevlist.push_back(*prev);
      while ( tmp != m_statusLine->begin() ) 
      {
	--tmp;
	if ( DoCurvesOverlap(*prev, *tmp))
	  prevlist.push_back(*tmp);
	else
	  break;
      }
 
      IntersectCurveGroup(*slIter, prevlist);
    }
    currentlist.push_back(*firstOne);

    EventCurveIter currentOne = firstOne; ++currentOne;
    EventCurveIter prevOne = firstOne;

    while ( currentOne != rightCurveEnd )
    {
      m_currentPos = m_sweepLinePos;
      PRINT_INSERT(*currentOne);
      ++slIter;
      slIter = m_statusLine->insert(slIter, *currentOne);
      SL_DEBUG(PrintStatusLine(););
      if ( DoCurvesOverlap(*currentOne, *prevOne))
      {
	IntersectCurveGroup(*currentOne, currentlist);
	currentlist.push_back(*currentOne);
      } else {
	prevlist = currentlist;
	currentlist.clear();
	currentlist.push_back(*currentOne);
      }
      
      IntersectCurveGroup(*currentOne, prevlist);
      prevOne = currentOne;
      ++currentOne;
    }

    lastOne = currentOne; --lastOne;
    m_currentPos = m_sweepLinePos;
    PRINT_INSERT(*lastOne);

    SL_DEBUG(PrintStatusLine();)
    StatusLineIter next = slIter; ++next;
    if ( next != m_statusLine->end() ) {
      IntersectCurveGroup(*next, currentlist);
      StatusLineIter tmp = next; ++tmp;
      while ( tmp != m_statusLine->end() ) 
      {
	if ( DoCurvesOverlap(*next, *tmp))
	{
	  IntersectCurveGroup(*tmp, currentlist);
	  ++tmp;
	}
	else
	  break;
      }
    }
  }
}



/*!
 * Perform intersection between the specified curve and all curves in the 
 * given group of curves.
 */ 
template <class CurveInputIterator, class SweepLineTraits_2,
         class SweepEvent, class CurveWrap>
inline void
Sweep_line_tight_2<CurveInputIterator, SweepLineTraits_2,SweepEvent,CurveWrap>::
IntersectCurveGroup(Subcurve *c1, SubCurveList &mylist)
{
  SL_DEBUG(std::cout << "Intersecting with " << mylist.size() << " curves\n";)
  SubCurveListIter i = mylist.begin();
  while ( i != mylist.end())
  {
    Intersect(c1, *i);
    ++i;
  }
}



/*! 
 * Finds intersection between two curves. 
 * If the two curves intersect, create a new event (or use the event that 
 * already exits in the intersection point) and insert the curves to the
 * event.
 * @param curve1 a pointer to the first curve
 * @param curve2 a pointer to the second curve
 * @return true if the two curves overlap.
*/
template <class CurveInputIterator, class SweepLineTraits_2,
         class SweepEvent, class CurveWrap>
inline bool 
Sweep_line_tight_2<CurveInputIterator, SweepLineTraits_2,SweepEvent,CurveWrap>::
Intersect(Subcurve *c1, Subcurve *c2)
{
  SL_DEBUG(std::cout << "Looking for intersection between:\n\t";)
  SL_DEBUG(c1->Print();)
  SL_DEBUG(std::cout << "\t";)
  SL_DEBUG(c2->Print();)
  SL_DEBUG(std::cout << "\n";)
	    
  if ( c1->getId() == c2->getId() ) {
    SL_DEBUG(std::cout << "same curve, returning....\n";)
    return false;
  }
  Subcurve *scv1 = c1;
  Subcurve *scv2 = c2;
  const X_curve_2 &cv1 = scv1->getCurve();
  const X_curve_2 &cv2 = scv2->getCurve();

  bool isOverlap = false;

  Point_2 p, p1;
  if ( m_traits->nearest_intersection_to_right(cv1, cv2, 
					       m_currentEvent->getPoint(), 
					       p, p1))
  {
    if ( !m_traits->point_is_same(p, p1)) {
      p = p1;
      SL_DEBUG(std::cout << "overlap detected\n";)
      isOverlap = true;
    }

    SL_DEBUG(
      std::cout << " a new event is created between:\n\t";
      scv1->Print();
      std::cout << "\t";
      scv2->Print();
      std::cout << "\trelative to ("
                << m_sweepLinePos << ")\n\t at (" 
                << p << ")" << std::endl;
    )

    // check to see if an event at this point already exists...
    EventQueueIter eqi = m_queue->find(p);
    Event *e = 0;
    if ( eqi == m_queue->end() )
    {
      e = new Event(p, m_traits); 
#ifndef NDEBUG
      e->id = m_eventId++;
#endif
      m_events.push_back(e);
      
      e->addCurveToLeft(c1, m_sweepLinePos);
      e->addCurveToLeft(c2, m_sweepLinePos);
      
      e->addCurveToRight(c1);
      e->addCurveToRight(c2);
      
      PRINT_NEW_EVENT(p, e);
      m_queue->insert(EventQueueValueType(p, e));
      return isOverlap;
    } else 
    {
      SL_DEBUG(std::cout << "event already exists, updating.. (" << p << ")\n";)
      e = eqi->second;
      if ( !scv1->isEndPoint(p)) {
	e->addCurveToLeft(c1, m_sweepLinePos);
	e->addCurveToRight(c1);
      }
      if ( !scv2->isEndPoint(p) ) {
	e->addCurveToLeft(c2, m_sweepLinePos);
	e->addCurveToRight(c2);
      }
      SL_DEBUG(e->Print();)
    }
    return isOverlap;
  } 
  SL_DEBUG(std::cout << "not found 2\n";)
  return isOverlap;
}



template <class CurveInputIterator,  class SweepLineTraits_2,
         class SweepEvent, class CurveWrap>
inline bool
Sweep_line_tight_2<CurveInputIterator, SweepLineTraits_2,SweepEvent,CurveWrap>::
isInternalXPoint(const Point_2 &p)
{
  EventListIter itt = m_miniq.begin();
  while ( itt != m_miniq.end() )
  {
    if ( m_traits->point_is_same(p, (*itt)->getPoint())) 
    {
      if ((*itt)->isInternalIntersectionPoint())
	return true;
      (*itt)->markInternalIntersectionPoint(); // this is to handle cases: |/ .
      return false;                            // (test 50/51)             |\ .
    } 
    ++itt;
  }
  assert(0);
  return false;
}



/*!
 * Handles the case in which a curve ont he status line passes through
 * one of the end points of the vertical curve.
 *
 * @param vcurve the vertical curve we are dealing with
 * @param curve a cerve that intersects with the vertical curve
 * @param topEndEvent the event attached to the top end of the vertical curve
 * @param tag 
 * @return returns true if the curve passed through one of the ends of the 
 *              vertical curve. Returns false otherwise.
 */template <class CurveInputIterator,  class SweepLineTraits_2,
         class SweepEvent, class CurveWrap>
inline bool
Sweep_line_tight_2<CurveInputIterator, SweepLineTraits_2,SweepEvent,CurveWrap>::
HandleVerticalCurveXAtEnd(Subcurve *vcurve, Subcurve *curve, 
			  Event *topEndEvent, SweepLineGetSubCurves tag)
{
  const Point_2 &topEnd = vcurve->getTopEnd();
  // handle a curve that goes through the top point of the vertical curve
  if (m_traits->curve_get_point_status(curve->getCurve(), topEnd) 
      == Traits::ON_CURVE )
  {
    if ( !curve->isLeftEnd(topEnd)) {
      topEndEvent->addCurveToLeft(curve, m_prevPos);
    }
    if ( ! curve->isRightEnd(topEnd)) {
      topEndEvent->addCurveToRight(curve);
    }
    return true;
  } 
  
  // handle a curve that goes through the bottom point of the vertical curve
  const Point_2 &currentPoint = m_currentEvent->getPoint();
  if (m_traits->curve_get_point_status((curve)->getCurve(), currentPoint) 
      == Traits::ON_CURVE)
  {
    if ( !(curve)->isLeftEnd(currentPoint)) {
      m_currentEvent->addCurveToLeft(curve, m_prevPos);
    }
    if ( ! (curve)->isRightEnd(currentPoint)) {
      m_currentEvent->addCurveToRight(curve);
    }
    return true;;
  }
  return false;
}



/*!
 * Handles the case in which a curve ont he status line passes through
 * one of the end points of the vertical curve.
 *
 * @param vcurve the vertical curve we are dealing with
 * @param curve a cerve that intersects with the vertical curve
 * @param topEndEvent the event attached to the top end of the vertical curve
 * @param tag 
 * @return returns true if the curve passed through one of the ends of the 
 *              vertical curve. Returns false otherwise.
 */
template <class CurveInputIterator,  class SweepLineTraits_2,
         class SweepEvent, class CurveWrap>
inline bool
Sweep_line_tight_2<CurveInputIterator, SweepLineTraits_2,SweepEvent,CurveWrap>::
HandleVerticalCurveXAtEnd(Subcurve *vcurve, Subcurve *curve, 
			  Event *topEndEvent, SweepLineGetPoints tag)
{
  const Point_2 &topEnd = vcurve->getTopEnd();
  // handle a curve that goes through the top point of the vertical curve
  if (m_traits->curve_get_point_status((curve)->getCurve(), topEnd) 
      == Traits::ON_CURVE )
  {
    if ( !curve->isEndPoint(topEnd))
      topEndEvent->markInternalIntersectionPoint();
    return true;
  } 

  // handle a curve that goes through the bottom point of the vertical curve
  if (m_traits->curve_get_point_status((curve)->getCurve(),
				       m_currentEvent->getPoint()) 
      == Traits::ON_CURVE)
  {
    if ( !curve->isEndPoint(m_currentEvent->getPoint()))
      m_currentEvent->markInternalIntersectionPoint();
    return true;
  }
  return false;
}


template <class CurveInputIterator,  class SweepLineTraits_2,
         class SweepEvent, class CurveWrap>
inline bool
Sweep_line_tight_2<CurveInputIterator, SweepLineTraits_2,SweepEvent,CurveWrap>::
DoCurvesOverlap(Subcurve *c1, Subcurve *c2)
{
  if ( m_traits->curve_compare_at_x_right(c1->getCurve(),
				    c2->getCurve(),
				    m_sweepLinePos) != EQUAL )
    return false;

  if ( m_traits->curves_overlap(c1->getCurve(),c2->getCurve()) )
    return true;

  return false;
}

template <class CurveInputIterator,  class SweepLineTraits_2,
         class SweepEvent, class CurveWrap>
inline bool
Sweep_line_tight_2<CurveInputIterator, SweepLineTraits_2,SweepEvent,CurveWrap>::
SimilarCurves(const X_curve_2 &a, const X_curve_2 &b)
{
  if ( m_traits->curve_is_same(a, b))
    return true;
  if ( m_traits->curve_is_same(m_traits->curve_flip(a), b))
    return true;
  return false;
}

template <class CurveInputIterator,  class SweepLineTraits_2,
          class SweepEvent, class CurveWrap>
inline bool
Sweep_line_tight_2<CurveInputIterator, SweepLineTraits_2,SweepEvent,CurveWrap>::
VerticalSubCurveExists(const X_curve_2 &a)
{
  for ( std::list<X_curve_2>::iterator iter = m_verticalSubCurves.begin() ;
	iter != m_verticalSubCurves.end() ; ++iter)
  {
    if (SimilarCurves(*iter, a)) 
      return true;
  }
  return false;
}



////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

template <class CurveInputIterator,  class SweepLineTraits_2,
         class SweepEvent, class CurveWrap>
inline void 
Sweep_line_tight_2<CurveInputIterator, SweepLineTraits_2,SweepEvent,CurveWrap>::
PrintEventQueue()
{
  SL_DEBUG(std::cout << std::endl << "Event queue: " << std::endl;)
  EventQueueIter iter = m_queue->begin();
  while ( iter != m_queue->end() )
  {
    SL_DEBUG(std::cout << "Point (" << iter->first << ")" << std::endl;)
    Event *e = iter->second;
    e->Print();
    ++iter;
  }
  SL_DEBUG(std::cout << "--------------------------------" << std::endl;)
}

template <class CurveInputIterator,  class SweepLineTraits_2,
          class SweepEvent, class CurveWrap>
inline void 
Sweep_line_tight_2<CurveInputIterator, SweepLineTraits_2,SweepEvent,CurveWrap>::
PrintSubCurves()
{
  SL_DEBUG(std::cout << std::endl << "Sub curves: " << std::endl;)
  SubCurveListIter iter = m_subCurves.begin();
  while ( iter != m_subCurves.end() )
  {
    (*iter)->Print();
    ++iter;
  }
}

template <class CurveInputIterator,  class SweepLineTraits_2,
         class SweepEvent, class CurveWrap>
inline void 
Sweep_line_tight_2<CurveInputIterator, SweepLineTraits_2,SweepEvent,CurveWrap>::
PrintStatusLine()
{
  if ( m_statusLine->size() == 0) {
    std::cout << std::endl << "Status line: empty" << std::endl;
    return;
  }
  std::cout << std::endl << "Status line: (" 
	    << m_currentPos << ")" << std::endl;
  StatusLineIter iter = m_statusLine->begin();
  while ( iter != m_statusLine->end() )
  {
    (*iter)->Print();
    ++iter;
  }
  std::cout << "Status line - end" << std::endl;
}

template <class CurveInputIterator,  class SweepLineTraits_2,
         class SweepEvent, class CurveWrap>
inline void 
Sweep_line_tight_2<CurveInputIterator, SweepLineTraits_2,SweepEvent,CurveWrap>::
PrintVerticals()
{
  if ( m_verticals.size() == 0) {
    std::cout << std::endl << "Verticals: empty" << std::endl;
    return;
  }
  std::cout << std::endl << "Verticals: " << m_verticals.size() << " (" 
	    << m_currentEvent->getPoint() << ")" << std::endl;
  SubCurveListIter iter = m_verticals.begin();
  while ( iter != m_verticals.end() )
  {
    (*iter)->Print();
    ++iter;
  }
  std::cout << "Verticals - end" << std::endl;
}

CGAL_END_NAMESPACE

#endif // CGAL_SWEEP_LINE_BASE_H
