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
// file          : include/CGAL/Sweep_line_base_2.h
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

#include <CGAL/Sweep_line_2/Sweep_line_event.h>
#include <CGAL/Sweep_line_2/Sweep_line_functors.h>
#include <CGAL/Sweep_line_2/Sweep_line_subcurve.h>

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
  Sweep_line_base_2 is a class that implements the sweep line algorithm
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



template <class CurveInputIterator,  class SweepLineTraits_2>
class Sweep_line_base_2
{
public:
  typedef SweepLineTraits_2 Traits;
  typedef typename Traits::Point_2 Point_2;
  typedef typename Traits::Curve_2 Curve_2;
  typedef typename Traits::X_curve_2 X_curve_2;

  typedef Sweep_line_event<Traits> Event;
  typedef Point_less_functor<Point_2, Traits> PointLess;
  typedef std::map<const Point_2, Event*, PointLess> EventQueue;
  typedef typename EventQueue::iterator EventQueueIter;
  typedef typename EventQueue::value_type EventQueueValueType;
  typedef std::vector<Event*> EventPtrContainer;
  typedef typename EventPtrContainer::iterator EventPtrContainerIter;

  typedef typename Event::SubCurveIter EventCurveIter;

  typedef Sweep_line_subcurve<Traits> SubCurve;
  typedef std::vector<SubCurve*> SubcurveContainer;
  typedef typename SubcurveContainer::iterator SubCurveIter;

  typedef Status_line_curve_less_functor<Traits> StatusLineCurveLess;
  typedef std::set<SubCurve*, StatusLineCurveLess> StatusLine;
  typedef typename StatusLine::iterator StatusLineIter;

  class  SweepLineGetSubCurves {};
  class  SweepLineGetPoints {};
  class  SweepLineGetInterCurveList {};

  Sweep_line_base_2() : m_traits(new Traits()), m_traitsOwner(true) {}
  Sweep_line_base_2(Traits *t) : m_traits(t), m_traitsOwner(false) {}

  ~Sweep_line_base_2();

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
  void  get_subcurves(CurveInputIterator begin, 
		      CurveInputIterator end, 
		      OutpoutIterator subcurves)
  { 
    Init(begin, end);
    SL_DEBUG(
      PrintSubCurves();
      PrintEventQueue();
    )
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
    Looping over the events in the queue, for each event we first
    handle the curves that are tothe left of the event point (i.e., 
    curves that we are done with), and then we look at the curves 
    to the right of the point, which means we attept to find intersections
    between them and their neighbours on the sweep line.
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
        HandleLeftCurves(out, tag);

        m_miniq.push_back(m_currentEvent);
        ++eventIter;
      }
      m_queue->erase(m_queue->begin(), eventIter);

      typename std::list<Event *>::iterator itt = m_miniq.begin();
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
  }

  /*! The main loop to calculate intersections among the curves
    Looping over the events in the queue, for each event we first
    handle the curves that are tothe left of the event point (i.e., 
    curves that we are done with), and then we look at the curves 
    to the right of the point, which means we attept to find intersections
    between them and their neighbours on the sweep line.
  */
  template <class OutpoutIterator>
  void Sweep(bool includeEndPoints, OutpoutIterator out,
	     SweepLineGetPoints tag)
  {
    //std::list<Event *> miniq;
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
	HandleLeftCurves(includeEndPoints, out, tag);

        m_miniq.push_back(m_currentEvent);
        ++eventIter;
      }

      m_queue->erase(m_queue->begin(), eventIter);

      typename std::list<Event *>::iterator itt = m_miniq.begin();
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


  /*!
    Handle a vertical curve when the event being processed is the top end 
    of the curve. In this situation, the event contains a list of intersection
    points on the vertical curve. We go through this list and outpt the 
    subcurves induced by these intersection points.
    If the curve is not vertical, returns without doing anything.

    @param out an iterator to the output
    @param tag a tag that indicates the version of the method
   */
  template <class OutpoutIterator>
  void HandleVerticalCurveTop(OutpoutIterator out, SweepLineGetSubCurves &tag)
  {
    SL_DEBUG(std::cout<<"HandleVerticalCurveTop... ";)
    if ( !m_currentEvent->doesContainVerticalCurve() ) {
      SL_DEBUG(std::cout<<"exiting\n ";)
      return;
    }
    SL_DEBUG(std::cout<<"\n ";)

    SubCurve *vcurve = m_currentEvent->getVerticalCurve();
    const Point_2 &topPoint = m_currentEvent->getPoint();
    // if this is the bottom point, nothing to do here
    if ( vcurve->isBottomEnd(topPoint)) {
      SL_DEBUG(std::cout<<"this is the bottom. exiting.\n";)
      return;
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
	      m_traits->curve_get_point_status((*slIter)->getCurve(), topPoint) 
	      == Traits::ABOVE_CURVE &&
	      m_traits->curve_get_point_status((*slIter)->getCurve(), 
					       vcurve->getBottomEnd()) 
	      == Traits::UNDER_CURVE )
      {
	SL_DEBUG(std::cout<<"checking \n";)
	SL_DEBUG((*slIter)->Print();) 
	if ( m_traits->compare_x((*slIter)->getLeftEnd(), topPoint) == EQUAL)
	{
	  m_currentEvent->addVerticalCurveXPoint((*slIter)->getLeftEnd());
	}
	++slIter;
      }   
      
    }

    // now we go over the list of intersection points on the vertical
    // curve in at the event and process them...
    SL_DEBUG(std::cout<<"handling the splitting now\n";)
    std::list<Point_2> &pointList = 
      m_currentEvent->getVerticalCurveXPointList();
    if ( pointList.empty() )
    {
      *out = m_currentEvent->getVerticalCurve()->getCurve();
      ++out;
      return;
    }
    
    X_curve_2 a, b, c;
    a = vcurve->getCurve();
    for ( std::list<Point_2>::iterator i = pointList.begin() ;
	  i != pointList.end(); ++i )
    {
      SL_DEBUG(std::cout<< "splitting: " << a << " at " << *i << "\n";)
      m_traits->curve_split(a, b, c, *i);
      if ( vcurve->isSourceLeftToTarget()) {
	*out = b; ++out;
	a = c;
      } else {
	*out = c; ++out;
	a = b;
      }
    }
    if ( vcurve->isSourceLeftToTarget() ) {
      *out = c; ++out;
    }
    else {
      *out = b; ++out;
    }
  }

  /*!
    Handle a vertical curve when the event being processed is the top end 
    of the curve. In this situation, the event contains a list of intersection
    points on the vertical curve. We go through this list and output the 
    intersection points.
    If the curve is not vertical, returns without doing anything.

    @param out an iterator to the output
    @param tag a tag that indicates the version of the method
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

    SubCurve *vcurve = m_currentEvent->getVerticalCurve();
    const Point_2 &topPoint = m_currentEvent->getPoint();
    if ( vcurve->isBottomEnd(topPoint)) {
      SL_DEBUG(std::cout<<"this is the bottom. exiting.\n";)
      return;
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
  }
 
  /*! For each left-curve, if it is the "last" subcurve, i.e., the 
    event point is the right-edge of the original curve, the 
    last sub curve is created and added to the result. Otherwise
    the curve is added as is to the result.
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
      SubCurve *leftCurve = *leftCurveIter; 
      const X_curve_2 &cv = leftCurve->getCurve();
      const Point_2 &lastPoint = leftCurve->getLastPoint();

      if ( leftCurve->isSource(eventPoint))
      {
        if ( !leftCurve->isTarget(lastPoint) )
        {
          X_curve_2 a,b;
          m_traits->curve_split(cv, a, b, lastPoint);
          *out = a; ++out;
        } else {
          *out = cv; ++out;
        }
      } else if ( leftCurve->isTarget(eventPoint))
      {
        if ( !leftCurve->isSource(lastPoint))
        {
          X_curve_2 a,b;
          m_traits->curve_split(cv, a, b, lastPoint);
          *out = b; ++out;
        } else {
          *out = cv; ++out;
        }

      } else { 
        X_curve_2 a,b;
        if ( leftCurve->isSource(lastPoint)) {
          m_traits->curve_split(cv, a, b, eventPoint);
          *out = a; ++out;
        } else if ( leftCurve->isTarget(lastPoint)) {
          m_traits->curve_split(cv, b, a, eventPoint);
          *out = a; ++out;
        } else {
          const X_curve_2 &lastCurve = leftCurve->getLastCurve();
          if ( leftCurve->isSourceLeftToTarget() ) {
            m_traits->curve_split(lastCurve, a, b, eventPoint);
            *out = a; ++out;
          } else {
            m_traits->curve_split(lastCurve, b, a, eventPoint);
            *out = a; ++out;
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
    event point is the right-edge of the original curve, the 
    last sub curve is created and added to the result. Otherwise
    the curve is added as is to the result.
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
      if ( includeEndPoints || m_currentEvent->isInternalIntersectionPoint())
        *out = eventPoint; ++out;    
      return;
    }
    EventCurveIter leftCurveIter = m_currentEvent->leftCurvesBegin();
    m_currentPos = m_prevPos;
    while ( leftCurveIter != m_currentEvent->leftCurvesEnd() )  // ** fix here
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
      *out = eventPoint; ++out;
  }

  

  void HandleRightCurves();
  Event *Intersect(SubCurve *c1, SubCurve *c2);
  Event *Intersect(SubCurve *c1, SubCurve *c2, SubCurve *c3);

  bool isInternalXPoint(const Point_2 &p);
  bool HandleVerticalCurveXAtEnd(SubCurve *vcurve, SubCurve *curve, 
				 Event *topEndEvent, SweepLineGetSubCurves tag);
  bool HandleVerticalCurveXAtEnd(SubCurve *vcurve, SubCurve *curve, 
				 Event *topEndEvent, SweepLineGetPoints tag);

  /*!
    When a curve is removed from the status line for good, its top and
    bottom neighbors become neighbors. This method finds these cases and
    looks for the itnersection point, if one exists.
    @param leftCurve a pointer to the curve that is about to be deleted
    @return an iterator to the position where the curve will be removed from.
  */
  StatusLineIter IntersectNeighboursAfterRemoval(SubCurve *leftCurve)
  {
    m_currentPos = m_prevPos;
    StatusLineIter sliter = m_statusLine->find(leftCurve);
    assert(sliter!=m_statusLine->end());
    StatusLineIter end = m_statusLine->end(); --end;
    if ( leftCurve->isEndPoint(m_currentEvent->getPoint()) && 
	 sliter != m_statusLine->begin() && sliter != end) 
    {
      StatusLineIter prev = sliter; --prev;
      StatusLineIter next = sliter; ++next;
      Intersect(*prev, *next);
    }
    return sliter;
  } 

  void PrintEventQueue();
  void PrintSubCurves();
  void PrintStatusLine();

private:
  /*! a pointer to a traits object */
  Traits *m_traits;

  /*! an indication to whether the traits should be deleted in the distructor */
  bool m_traitsOwner;

  /*! used to hold all event, processed or not, so that they can 
`     all be deleted at the end */
  EventPtrContainer m_events; 

  /*! the queue of events (intersection points) to handle */
  EventQueue *m_queue;

  /*! The subcurves, as created on the fly */
  SubcurveContainer m_subCurves;

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
  typename std::list<Event *> m_miniq;

#ifndef NDEBUG
  int m_eventId;
  int m_curveId;
#endif
};

template <class CurveInputIterator,  class SweepLineTraits_2>
Sweep_line_base_2<CurveInputIterator, SweepLineTraits_2>::
~Sweep_line_base_2() 
{
  if ( m_traitsOwner ) delete m_traits;

  for ( SubCurveIter sci = m_subCurves.begin() ; 
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
  - x-monotonize the inf\put curves
  - for each end point of each curve create an event
  - initialize the event queue
  -
*/
template <class CurveInputIterator,  class SweepLineTraits_2>
inline void 
Sweep_line_base_2<CurveInputIterator, SweepLineTraits_2>::
Init(CurveInputIterator begin, CurveInputIterator end)
{
  PointLess pred(m_traits);
  m_queue = new EventQueue(pred);
  StatusLineCurveLess slcurveless(m_traits);
  m_statusLine = new StatusLine(slcurveless);

  SL_DEBUG(m_eventId = 0;)
  SL_DEBUG(m_curveId = 0 ;)

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
    SL_DEBUG(m_curveId++;)
  }
}

/*! Given an x-monotone curve, create events for each end (if 
    one doesn't exist already). 
    For each curve create a SubCurve instance.
 */
template <class CurveInputIterator,  class SweepLineTraits_2>
inline void 
Sweep_line_base_2<CurveInputIterator, SweepLineTraits_2>::
InitCurve(X_curve_2 &curve)
{
  const Point_2 &source = m_traits->curve_source(curve);
  const Point_2 &target = m_traits->curve_target(curve);
  Event *e = 0;
  
  SubCurve *subCv = new SubCurve(curve, &m_currentPos, m_traits);
  SL_DEBUG(subCv->id = m_curveId;)
  m_subCurves.push_back(subCv);
  
  // handle the source point
  EventQueueIter eventIter = m_queue->find(source);
  if ( eventIter != m_queue->end() ) {
    SL_DEBUG(std::cout << "event " << source << " already exists\n";)
    e = eventIter->second;
    e->markInternalIntersectionPoint();
  } else  {
    e = new Event(source, m_traits); 
    SL_DEBUG(e->id = m_eventId++;)
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
    e->markInternalIntersectionPoint();
  } else  {
    e = new Event(target, m_traits); 
    SL_DEBUG(e->id = m_eventId++;)
    m_events.push_back(e);
    m_queue->insert(EventQueueValueType(target, e));
  }
  e->addCurve(subCv);
  PRINT_NEW_EVENT(target, e);
}

/*! This pass comes to take care of cases in which we reach an
    event point that is the left end of a curve that starts on another curve.
    For example:
         /      ------         /
        /----      \          /
       /            \     ----------
    This method is called before the left curves of an event are handled.
*/
template <class CurveInputIterator,  class SweepLineTraits_2>
inline void 
Sweep_line_base_2<CurveInputIterator, SweepLineTraits_2>::
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
 
    if ( slIter != m_statusLine->begin() )
    { 
      --prev;
      if ( m_traits->curve_get_point_status((*prev)->getCurve(), p) ==
	   Traits::ON_CURVE && !(*prev)->isEndPoint(p))
      {
	std::cout << "here\n";
	m_currentEvent->addCurveToRight(*prev);
	m_currentEvent->addCurveToLeft(*prev, m_prevPos);
	m_currentEvent->markInternalIntersectionPoint();
      }
    }
    if ( slIter != m_statusLine->end() )
    {
      if ( m_traits->curve_get_point_status((*next)->getCurve(), p) ==
	   Traits::ON_CURVE && !(*next)->isEndPoint(p))
      {
	std::cout << "here\n";
	m_currentEvent->addCurveToRight(*next);
	m_currentEvent->addCurveToLeft(*next, m_prevPos);
	m_currentEvent->markInternalIntersectionPoint();
      }    
    } 
    ++rightIter;
  }

  SL_DEBUG(std::cout << "First pass - done\n" ;)
}

/*!
  Handles the degenerate case of vertical curves. Most of the cases
  that occur with vertical curves are handled by this method and 
  HandleVerticalCurveTop method.
  When the current event is the bottom end of a vertical curve, we look
  for intersection points between the vertical curve and any curve
  in the status line that in the y-range that is defined by the bottom 
  and top ends of the vertical curve. When those are found, we create
  new events, unless ones already exist, in which case we update the events.
  
  @param tag a tag that indicates the version of this method
  \sa HandleVerticalCurveTop
 */
template <class CurveInputIterator,  class SweepLineTraits_2>
inline void 
Sweep_line_base_2<CurveInputIterator, SweepLineTraits_2>::
HandleVerticalCurveBottom(SweepLineGetSubCurves &tag)
{
  SL_DEBUG(std::cout<<"HandleVerticalCurveBottom... ";)
  if ( !m_currentEvent->doesContainVerticalCurve() )
  {
    SL_DEBUG(std::cout<<" - not vertical - exiting\n ";)
    return;
  }
  SL_DEBUG(std::cout<<"\n ";)

  SubCurve *vcurve = m_currentEvent->getVerticalCurve();
  const Point_2 &currentPoint = m_currentEvent->getPoint();
  if ( vcurve->isTopEnd(currentPoint))
    return;

  SL_DEBUG(std::cout<<"handling bottom point of vertical curve\n";)
  StatusLineIter slIter = m_statusLine->lower_bound(vcurve);
  if ( slIter == m_statusLine->end() ) {
    SL_DEBUG(std::cout<<"no curves intersecting. exiting\n";)
    return;
  }    

  SL_DEBUG(std::cout<<"starting at curve \n";)
  SL_DEBUG((*slIter)->Print();)
  const Point_2 &topEnd = vcurve->getTopEnd();
  EventQueueIter topEndEventIter = m_queue->find(topEnd);
  assert(topEndEventIter!=m_queue->end());
  Event *topEndEvent = topEndEventIter->second;

  std::cout << "before while\n";
  while ( slIter != m_statusLine->end() &&
	  m_traits->curve_get_point_status((*slIter)->getCurve(), topEnd) 
	  != Traits::UNDER_CURVE &&
	  m_traits->curve_get_point_status((*slIter)->getCurve(), 
					   currentPoint) 
	  != Traits::ABOVE_CURVE )
  {
    SL_DEBUG(std::cout<<"intersecting with \n";)
    SL_DEBUG((*slIter)->Print();) 
    const Point_2 &currentPoint = m_currentEvent->getPoint();

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
    assert(res==true);
    EventQueueIter eqi = m_queue->find(p);
    Event *e = 0;
    if ( eqi == m_queue->end() )
    {
      e = new Event(p, m_traits, true); 
      SL_DEBUG(e->id = m_eventId++;)
      m_events.push_back(e);
      
      e->addCurveToLeft(*slIter, m_sweepLinePos);
      e->addCurveToRight(*slIter);
      PRINT_NEW_EVENT(p, e);
      m_queue->insert(EventQueueValueType(p, e));

    } else {

      e = eqi->second;
      e->markInternalIntersectionPoint();
      SL_DEBUG(std::cout << "Updating event \n";)
      SL_DEBUG(e->Print();)
    }

    topEndEvent->addVerticalCurveXPoint(p);
    ++slIter;
  }    
  std::cout << "after while \n";

  SL_DEBUG(std::cout<<"Done Handling vertical\n";)
}

/*!
  Handles the degenerate case of vertical curves. Most of the cases
  that occur with vertical curves are handled by this method and the
  HandleVerticalCurveTop method.

  When the current event is the bottom end of a vertical curve, we look
  for intersection points between the vertical curve and any curve
  in the status line that in the y-range that is defined by the bottom 
  and top ends of the vertical curve. When those are found, we create
  new events, unless ones already exist, in which case we update the events.
  
  @param tag a tag that indicates the version of this method
  \sa HandleVerticalCurveTop
 */
template <class CurveInputIterator,  class SweepLineTraits_2>
inline void 
Sweep_line_base_2<CurveInputIterator, SweepLineTraits_2>::
HandleVerticalCurveBottom(SweepLineGetPoints &tag)
{
  SL_DEBUG(std::cout<<"HandleVerticalCurveBottom... ";)
  if ( !m_currentEvent->doesContainVerticalCurve() )
  {
    SL_DEBUG(std::cout<<"exiting\n ";)
    return;
  }
  SL_DEBUG(std::cout<<"\n ";)

  SubCurve *vcurve = m_currentEvent->getVerticalCurve();
  const Point_2 &currentPoint = m_currentEvent->getPoint();

  if ( vcurve->isTopEnd(currentPoint))
    return;

  SL_DEBUG(std::cout<<"handling bottom point of vertical curve\n";)
  StatusLineIter slIter = m_statusLine->lower_bound(vcurve);
  if ( slIter == m_statusLine->end() ) {
    SL_DEBUG(std::cout<<"no curves intersecting. exiting\n";)
    return;
  }    

  SL_DEBUG(std::cout<<"starting at curve \n";)
  SL_DEBUG((*slIter)->Print();)

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
    bool res = m_traits->nearest_intersection_to_right(vcurve->getCurve(), 
						       (*slIter)->getCurve(), 
						       currentPoint, 
						       xp, xp);
    assert(res==true);
    EventQueueIter eqi = m_queue->find(xp);
    Event *e = 0;
    if ( eqi == m_queue->end() )
    {
      e = new Event(xp, m_traits, true); 
      SL_DEBUG(e->id = m_eventId++;)
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

  SL_DEBUG(std::cout<<"Done Handling vertical\n";)
}


/*! Loop over the curves to the right of the sweep line and handle them:
  - if we are at the beginning of the curve, we insert it to the sweep 
    line, then we look if it intersects any of its neighbours.
  - if we are at an intersection point between two curves, we add them
    to the sweep line and attempt to intersect them with their neighbours. 
    We also check to see if the two intersect again to the right of the point.
 */
template <class CurveInputIterator,  class SweepLineTraits_2>
inline void 
Sweep_line_base_2<CurveInputIterator, SweepLineTraits_2>::
HandleRightCurves()
{
  SL_DEBUG(std::cout << "Handling right curves\n" ;)
  SL_DEBUG(std::cout << m_currentEvent->getPoint() << "\n";)
  int numRightCurves = m_currentEvent->getNumRightCurves();
  if ( numRightCurves == 0 )
    return;

  m_currentPos = m_sweepLinePos;
  if ( numRightCurves == 1 )
  {
    SL_DEBUG(std::cout << " - beginning of curve " << std::endl;)

    SL_DEBUG(
      SubCurve *tmp1 = *(m_currentEvent->rightCurvesBegin());
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
    next++;

    if ( next == m_statusLine->end() )
    { 
      --prev;
      Intersect(*slIter, *prev);
    } else if ( slIter == m_statusLine->begin() )
    {
       Intersect(*next, *slIter);
    } else {

      --prev;
      SubCurve *prevSC = *prev;
      SubCurve *SC = *slIter;
      SubCurve *nextSC = *next;
      Intersect(prevSC, SC, nextSC);
    }
 
  } else
  {

    SL_DEBUG(std::cout << " - intersection point " << std::endl;)
    EventCurveIter firstOne = m_currentEvent->rightCurvesBegin();
    EventCurveIter lastOne = m_currentEvent->rightCurvesEnd(); --lastOne;

    PRINT_INSERT(*firstOne);
    std::pair<StatusLineIter, bool> tmp = m_statusLine->insert(*firstOne);
    StatusLineIter slIter = tmp.first;

    SL_DEBUG(PrintStatusLine();)
    StatusLineIter prev = slIter;
    if ( slIter != m_statusLine->begin() )
    { 
      --prev;
      Intersect(*prev, *slIter);
    }

    prev = slIter;
    EventCurveIter currentOne = firstOne; ++currentOne;

    while ( currentOne != lastOne )
    {
      m_currentPos = m_sweepLinePos;
      PRINT_INSERT(*currentOne);
      ++slIter;
      slIter = m_statusLine->insert(slIter, *currentOne);
      SL_DEBUG(PrintStatusLine();)
      Intersect(*prev, *slIter);
      prev = slIter;
      ++currentOne;
    }

    Intersect(*prev, *lastOne);

    m_currentPos = m_sweepLinePos;
    PRINT_INSERT(*lastOne);
    ++slIter;
    slIter = m_statusLine->insert(slIter, *lastOne);
    SL_DEBUG(PrintStatusLine();)
    StatusLineIter next = slIter; ++next;
    if ( next != m_statusLine->end() ) {
      Intersect(*slIter, *next);
    }
  }
}


/*! 
  Finds intersection between two curves. 
  @param curve1 a pointer to the first curve
  @param curve2 a pointer to the second curve
  @return a pointer to the event. 0 if it does not exist.
*/
template <class CurveInputIterator, class SweepLineTraits_2>
inline Sweep_line_event<SweepLineTraits_2> *
Sweep_line_base_2<CurveInputIterator, SweepLineTraits_2>::
Intersect(SubCurve *c1, SubCurve *c2)
{
  SL_DEBUG(std::cout << "Looking for intersection between:\n\t";)
  SL_DEBUG(c1->Print();)
  SL_DEBUG(std::cout << "\t";)
  SL_DEBUG(c2->Print();)
  SL_DEBUG(std::cout << "\n";)
	    
  SubCurve *scv1 = c1;
  SubCurve *scv2 = c2;
  const X_curve_2 &cv1 = scv1->getCurve();
  const X_curve_2 &cv2 = scv2->getCurve();

  Point_2 p;
  if ( m_traits->nearest_intersection_to_right(cv1, cv2, 
					       m_currentEvent->getPoint(), 
					       p, p))
  {
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
      e = new Event(p, m_traits, true); 
      SL_DEBUG(e->id = m_eventId++;)
      m_events.push_back(e);
      
      //m_currentPos = m_sweepLinePos;
      e->addCurveToLeft(c1, m_sweepLinePos);
      e->addCurveToLeft(c2, m_sweepLinePos);
      
      e->addCurveToRight(c1);
      e->addCurveToRight(c2);
      
      PRINT_NEW_EVENT(p, e);
      m_queue->insert(EventQueueValueType(p, e));
      return e;
    } else 
    {
      SL_DEBUG(std::cout << "event already exists, updating.. (" << p << ")\n";)
      e = eqi->second;
      e->markInternalIntersectionPoint();
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
    return e;
  } 
  SL_DEBUG(std::cout << "not found 2\n";)

  return 0;
  
}

template <class CurveInputIterator,  class SweepLineTraits_2>
inline Sweep_line_event<SweepLineTraits_2> *
Sweep_line_base_2<CurveInputIterator, SweepLineTraits_2>::
Intersect(SubCurve *c1, SubCurve *c2, SubCurve *c3)
{
  const X_curve_2 &cv1 = c1->getCurve();
  const X_curve_2 &cv2 = c2->getCurve();
  const X_curve_2 &cv3 = c3->getCurve();
  Point_2 p12, p23;

  bool _12intersect = 
    m_traits->nearest_intersection_to_right(cv1, cv2, m_sweepLinePos, p12, p12);
  bool _23intersect = 
    m_traits->nearest_intersection_to_right(cv2, cv3, m_sweepLinePos, p23, p23);

  if ( !_12intersect && ! _23intersect ){
    return 0;
  }

  Point_2 p;
  SubCurve *a, *b;
  if ( _12intersect && ! _23intersect )
  {
    a = c1; b = c2;
    p = p12;

  } else if ( !_12intersect && _23intersect )
  {
    a = c2; b = c3;
    p = p23;

  } else {

    if ( m_traits->compare_x(p12, p23) == SMALLER )
    {
      p = p12;
      a = c1; b = c2;
    } else {
      p = p23;
      a = c2; b = c3;
    }
  }

  // check to see if an event at this point already exists...
  EventQueueIter eqi = m_queue->find(p);
  Event *e = 0;
  if ( eqi == m_queue->end() )
  {
    e = new Event(p, m_traits, true); 
    SL_DEBUG(e->id = m_eventId++;)
    m_events.push_back(e);
    
    e->addCurveToLeft(a, m_sweepLinePos);
    e->addCurveToLeft(b, m_sweepLinePos);
    
    e->addCurveToRight(a);
    e->addCurveToRight(b);
    
    PRINT_NEW_EVENT(p, e);
    m_queue->insert(EventQueueValueType(p, e));

  } else {

    e = eqi->second;
    e->markInternalIntersectionPoint();
    if ( !a->isEndPoint(p))
    {
      e->addCurveToLeft(a, m_sweepLinePos);
      e->addCurveToRight(a);
    }
    if ( !b->isEndPoint(p) )
    {
      e->addCurveToLeft(b, m_sweepLinePos);
      e->addCurveToRight(b);
    }
    SL_DEBUG(e->Print();)
  }
  return e;
  
}




template <class CurveInputIterator,  class SweepLineTraits_2>
inline bool
Sweep_line_base_2<CurveInputIterator, SweepLineTraits_2>::
isInternalXPoint(const Point_2 &p)
{
  typename std::list<Event *>::iterator itt = m_miniq.begin();
  while ( itt != m_miniq.end() )
  {
    if ( m_traits->point_is_same(p, (*itt)->getPoint())) 
    {
      if ((*itt)->isInternalIntersectionPoint())
	return true;
      return false;
    } 
    ++itt;
  }
  assert(0);
}

/*!
  Handles the case in which a curve ont he status line passes through
  one of the end points of the vertical curve.

  @param vcurve the vertical curve we are dealing with
  @param curve a cerve that intersects with the vertical curve
  @param topEndEvent the event attached to the top end of the vertical curve
  @param tag 
  @return returns true if the curve passed through one of the ends of the 
               vertical curve. Returns false otherwise.
 */template <class CurveInputIterator,  class SweepLineTraits_2>
inline bool
Sweep_line_base_2<CurveInputIterator, SweepLineTraits_2>::
HandleVerticalCurveXAtEnd(SubCurve *vcurve, SubCurve *curve, 
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
  Handles the case in which a curve ont he status line passes through
  one of the end points of the vertical curve.

  @param vcurve the vertical curve we are dealing with
  @param curve a cerve that intersects with the vertical curve
  @param topEndEvent the event attached to the top end of the vertical curve
  @param tag 
  @return returns true if the curve passed through one of the ends of the 
               vertical curve. Returns false otherwise.
 */
template <class CurveInputIterator,  class SweepLineTraits_2>
inline bool
Sweep_line_base_2<CurveInputIterator, SweepLineTraits_2>::
HandleVerticalCurveXAtEnd(SubCurve *vcurve, SubCurve *curve, 
			  Event *topEndEvent, SweepLineGetPoints tag)
{
  const Point_2 &topEnd = vcurve->getTopEnd();
  // handle a curve that goes through the top point of the vertical curve
  if (m_traits->curve_get_point_status((curve)->getCurve(), topEnd) 
      == Traits::ON_CURVE )
  {
    topEndEvent->markInternalIntersectionPoint();
    return true;
  } 

  // handle a curve that goes through the bottom point of the vertical curve
  if (m_traits->curve_get_point_status((curve)->getCurve(),
				       m_currentEvent->getPoint()) 
      == Traits::ON_CURVE)
  {
    m_currentEvent->markInternalIntersectionPoint();
    return true;
  }
  return false;
}

////////////////////////////////////////////////////////////////////////////

template <class CurveInputIterator,  class SweepLineTraits_2>
inline void 
Sweep_line_base_2<CurveInputIterator, SweepLineTraits_2>::
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

template <class CurveInputIterator,  class SweepLineTraits_2>
inline void 
Sweep_line_base_2<CurveInputIterator, SweepLineTraits_2>::
PrintSubCurves()
{
  SL_DEBUG(std::cout << std::endl << "Sub curves: " << std::endl;)
  SubCurveIter iter = m_subCurves.begin();
  while ( iter != m_subCurves.end() )
  {
    (*iter)->Print();
    ++iter;
  }
}

template <class CurveInputIterator,  class SweepLineTraits_2>
inline void 
Sweep_line_base_2<CurveInputIterator, SweepLineTraits_2>::
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

CGAL_END_NAMESPACE

#endif // CGAL_SL_H
