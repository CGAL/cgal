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
// file          : include/CGAL/SL.h
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
    PerformIntersection(subcurves);

    SL_DEBUG(PrintResult();)
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
    PerformIntersection_for_points(includeEndPoints, points);
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
  }
protected:
  void Init(CurveInputIterator begin, CurveInputIterator end);

  void InitCurve(X_curve_2 &curve);
  void UpdateEventCurves(Event *e, SubCurve *subCurve);

  /*! The main loop to calculate intersections among the curves
    Looping over the events in the queue, for each event we first
    handle the curves that are tothe left of the event point (i.e., 
    curves that we are done with), and then we look at the curves 
    to the right of the point, which means we attept to find intersections
    between them and their neighbours on the sweep line.
  */
  template <class OutpoutIterator>
  void PerformIntersection(OutpoutIterator out)
  {
    std::list<Event *> miniq;
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
        SL_DEBUG(std::cout << "--------------------------------------"
                 << std::endl;
                 std::cout << "------------- " << *p << " --------------"
                 << std::endl;
                 std::cout << "--------------------------------------"
                 << std::endl;
                 PrintStatusLine();
                 m_currentEvent->Print();
                 std::cout << "--------------------------------------"
                 << std::endl;
                 std::cout << "--------------------------------------"
                 << std::endl;
                 )

        FirstPass();
        HandleLeftCurves_for_subcurves(*p, out);

        miniq.push_back(m_currentEvent);
        ++eventIter;
      }

      m_queue->erase(m_queue->begin(), eventIter);

      typename std::list<Event *>::iterator itt = miniq.begin();
      while ( itt != miniq.end())
      {
        m_currentEvent = *itt;
        HandleRightCurves();
        ++itt;
      }
      miniq.clear();
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
  void PerformIntersection_for_points(bool includeEndPoints,
					OutpoutIterator out)
  {
    std::list<Event *> miniq;
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
        SL_DEBUG(std::cout << "--------------------------------------"
                 << std::endl;
                 std::cout << "------------- " << *p << " --------------"
                 << std::endl;
                 std::cout << "--------------------------------------"
                 << std::endl;
                 PrintStatusLine();
                 m_currentEvent->Print();
                 std::cout << "--------------------------------------"
                 << std::endl;
                 std::cout << "--------------------------------------"
                 << std::endl;
        )

        FirstPass();
        HandleLeftCurves_for_points(*p, includeEndPoints, out);

        miniq.push_back(m_currentEvent);
        ++eventIter;
      }

      m_queue->erase(m_queue->begin(), eventIter);

      typename std::list<Event *>::iterator itt = miniq.begin();
      while ( itt != miniq.end())
      {
        m_currentEvent = *itt;
        HandleRightCurves();
        ++itt;
      }
      miniq.clear();
      eventIter = m_queue->begin();
    }
  }

  void FirstPass();

  /*! For each left-curve, if it is the "last" subcurve, i.e., the 
    event point is the right-edge of the original curve, the 
    last sub curve is created and added to the result. Otherwise
    the curve is added as is to the result.
  */
  template <class OutpoutIterator>
  void HandleLeftCurves_for_subcurves(const Point_2 &p, OutpoutIterator out)
  {
    SL_DEBUG(std::cout << "Handling left curve" << std::endl;)
    SL_DEBUG(m_currentEvent->Print();)
    EventCurveIter leftCurveIter = m_currentEvent->leftCurvesBegin();
    m_currentPos = m_prevPos;
    while ( leftCurveIter != m_currentEvent->leftCurvesEnd() )  // ** fix here
    {
      SubCurve *leftCurve = *leftCurveIter; 
      const X_curve_2 &cv = leftCurve->getCurve();
      const Point_2 &lastPoint = leftCurve->getLastPoint();

      if ( leftCurve->isSource(p))
      {
        if ( !leftCurve->isTarget(lastPoint) )
        {
          X_curve_2 a,b;
          m_traits->curve_split(cv, a, b, lastPoint);
          *out = a; ++out;
        } else {
          *out = cv; ++out;
        }
      } else if ( leftCurve->isTarget(p))
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
          m_traits->curve_split(cv, a, b, p);
          *out = a; ++out;
        } else if ( leftCurve->isTarget(lastPoint)) {
          m_traits->curve_split(cv, b, a, p);
          *out = a; ++out;
        } else {
          const X_curve_2 &lastCurve = leftCurve->getLastCurve();
          if ( leftCurve->isSourceLeftToTarget() ) {
            m_traits->curve_split(lastCurve, a, b, p);
            *out = a; ++out;
          } else {
            m_traits->curve_split(lastCurve, b, a, p);
            *out = a; ++out;
          }
        }
        leftCurve->setLastPoint(p);
        leftCurve->setLastCurve(b); 
      }
      m_currentPos = m_prevPos;
      SL_DEBUG(PrintStatusLine();)
      PRINT_ERASE(leftCurve);

      // when a curve is deleted from the sweep line for good (i.e., 
      // the sweepline position is at the end of the curve), and it 
      // is not the first or last on the sweep line, its top neighbour 
      // and bottom neighbour become neighbours and we need to check 
      // for intersection.
      StatusLineIter sliter = m_statusLine->find(leftCurve);
      assert(sliter!=m_statusLine->end());
      StatusLineIter end = m_statusLine->end(); --end;
      if ( leftCurve->isEndPoint(p) && 
           sliter != m_statusLine->begin() && sliter != end) 
      {
        StatusLineIter prev = sliter; --prev;
        StatusLineIter next = sliter; ++next;
        Intersect(*prev, *next);
      } 

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
  void HandleLeftCurves_for_points(const Point_2 &p, bool includeEndPoints,
                                   OutpoutIterator out)
  {
    SL_DEBUG(std::cout << "Handling left curve" << std::endl;)
    SL_DEBUG(m_currentEvent->Print();)
    if ( !m_currentEvent->hasLeftCurves() )
    {
      if ( includeEndPoints )
        *out = p; ++out;    
      return;
    }
    EventCurveIter leftCurveIter = m_currentEvent->leftCurvesBegin();
    m_currentPos = m_prevPos;
    bool isEndPoint = true;
    while ( leftCurveIter != m_currentEvent->leftCurvesEnd() )  // ** fix here
    {
      m_currentPos = m_prevPos;
      SL_DEBUG(PrintStatusLine();)
      SubCurve *leftCurve = *leftCurveIter; 
      isEndPoint &= leftCurve->isEndPoint(p);
      PRINT_ERASE(leftCurve);

      // when a curve is deleted from the sweep line for good (i.e., 
      // the sweepline position is at the end of the curve), and it 
      // is not the first or last on the sweep line, its top neighbour 
      // and bottom neighbour become neighbours and we need to check 
      // for intersection.
      StatusLineIter sliter = m_statusLine->find(leftCurve);
      assert(sliter!=m_statusLine->end());
      StatusLineIter end = m_statusLine->end(); --end;
      if ( leftCurve->isEndPoint(p) && 
           sliter != m_statusLine->begin() && sliter != end) 
      {
        StatusLineIter prev = sliter; --prev;
        StatusLineIter next = sliter; ++next;
        Intersect(*prev, *next);
      } 

      m_currentPos = m_prevPos;
      m_statusLine->erase(sliter);
      ++leftCurveIter;
    }

    if ( includeEndPoints || !isEndPoint )
      *out = p; ++out;
  }
  
  void HandleRightCurves();
  Event *Intersect(SubCurve *c1, SubCurve *c2);
  Event *Intersect(SubCurve *c1, SubCurve *c2, SubCurve *c3);

  void PrintEventQueue();
  void PrintSubCurves();
  void PrintStatusLine();
  void PrintResult();

private:
  Traits *m_traits;
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

  /*! A reference point that is used for comapring curves */
  Point_2 m_currentPos;

  /*! A reference point that is used for comapring curves */
  Point_2 m_prevPos;

  /*! The current position (in X) of the sweep line */
  Point_2 m_sweepLinePos;

  /*! the final subcurves */
  std::vector<X_curve_2> m_result;

  /*! if non x-monotone are specified, this hold the x-monotone 
    curves created from them */
  std::vector<X_curve_2> m_xcurves;

  Event *m_currentEvent;

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
  } else  {
    e = new Event(source, m_traits); 
    SL_DEBUG(e->id = m_eventId++;)
    m_events.push_back(e);
    m_queue->insert(EventQueueValueType(source, e));
  }
  UpdateEventCurves(e, subCv);
  PRINT_NEW_EVENT(source, e);
    
  // handle the target point
  eventIter = m_queue->find(target);
  if ( eventIter != m_queue->end() ) {
    SL_DEBUG(std::cout << "event " << target << " already exists\n";)
    e = eventIter->second;
  } else  {
    e = new Event(target, m_traits); 
    SL_DEBUG(e->id = m_eventId++;)
    m_events.push_back(e);
    m_queue->insert(EventQueueValueType(target, e));
  }
  UpdateEventCurves(e, subCv);
  PRINT_NEW_EVENT(target, e);
}

/*! Given an event (the point indicates the position of the event) 
    and a curve, this method adds
    the curve to the event according to its position relative to 
    the point (right or left, or both).
*/
template <class CurveInputIterator,  class SweepLineTraits_2>
inline void 
Sweep_line_base_2<CurveInputIterator, SweepLineTraits_2>::
UpdateEventCurves(Event *e, SubCurve *subCurve)
{
  const Point_2 &point = e->getPoint();
  const X_curve_2 &curve = subCurve->getCurve();
  const Point_2 &source = m_traits->curve_source(curve);
  const Point_2 &target = m_traits->curve_target(curve);

  if ( point == source || point == target ) 
  {
    const Point_2 *rel = &(source);
    if ( point == source )
      rel = &(target);
    
    
    if ( m_traits->compare_x(point, *rel) == LARGER ) {
      m_currentPos = *rel;
      e->addCurveToLeft(subCurve);
    } else {
      m_currentPos = point;
      e->addCurveToRight(subCurve);
    }

  } else {

    if ( m_sweepLinePos != point )
      m_currentPos = m_sweepLinePos;
    else 
      m_currentPos = m_prevPos;
    e->addCurveToLeft(subCurve);
    m_currentPos = point;
    e->addCurveToRight(subCurve);
  }
}

/*! This pass comes to take care of cases in which we reach an
    event point that is an end point of a curve, and also ends 
    on another curve:
              /
        -----/
            /
*/
template <class CurveInputIterator,  class SweepLineTraits_2>
inline void 
Sweep_line_base_2<CurveInputIterator, SweepLineTraits_2>::
FirstPass()
{
  if ( m_statusLine->size() == 0 )
    return;
  const Point_2& p = m_currentEvent->getPoint();

  DBG("First pass - right curves");
  EventCurveIter rightIter = m_currentEvent->rightCurvesBegin();
  m_currentPos = m_sweepLinePos;
  while ( rightIter != m_currentEvent->rightCurvesEnd())
  {
    // if the point is not an edge point, return
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
	m_currentPos = m_sweepLinePos;
	m_currentEvent->addCurveToRight(*prev);
	m_currentPos = m_prevPos;
	m_currentEvent->addCurveToLeft(*prev);
      }
    }
    if ( slIter != m_statusLine->end() )
    {
      if ( m_traits->curve_get_point_status((*next)->getCurve(), p) ==
	   Traits::ON_CURVE && !(*next)->isEndPoint(p))
      {
	m_currentPos = m_sweepLinePos;
	m_currentEvent->addCurveToRight(*next);
	m_currentPos = m_prevPos;
	m_currentEvent->addCurveToLeft(*next);
      }    
    } 
    ++rightIter;
  }

  SL_DEBUG(std::cout << "First pass - done\n" ;)
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
  int numRightCurves = m_currentEvent->getNumRightCurves();
  if ( numRightCurves == 0 )
    return;

  m_currentPos = m_sweepLinePos;
  if ( numRightCurves == 1 )
  {
    // this is the beginning of the curve, insert into the status line
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

// c1 is in the event already and c2 is not
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
      e = new Event(p, m_traits); 
      SL_DEBUG(e->id = m_eventId++;)
      m_events.push_back(e);
      
      m_currentPos = m_sweepLinePos;
      e->addCurveToLeft(c1);
      e->addCurveToLeft(c2);
      
      m_currentPos = p;
      e->addCurveToRight(c1);
      e->addCurveToRight(c2);
      
      PRINT_NEW_EVENT(p, e);
      m_queue->insert(EventQueueValueType(p, e));
      return e;
    } else 
    {
      SL_DEBUG(std::cout << "event already exists, updating.. (" << p << ")\n";)
      e = eqi->second;

      if ( !scv1->isEndPoint(p)) {
	m_currentPos = m_sweepLinePos;
	e->addCurveToLeft(c1);
	m_currentPos = p;
	e->addCurveToRight(c1);
      }
      if ( !scv2->isEndPoint(p) ) {
	m_currentPos = m_sweepLinePos;
	e->addCurveToLeft(c2);
	m_currentPos = p;
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
    e = new Event(p, m_traits); 
    SL_DEBUG(e->id = m_eventId++;)
    m_events.push_back(e);
    
    m_currentPos = m_sweepLinePos;
    e->addCurveToLeft(a);
    e->addCurveToLeft(b);
    
    m_currentPos = p;
    e->addCurveToRight(a);
    e->addCurveToRight(b);
    
    PRINT_NEW_EVENT(p, e);
    m_queue->insert(EventQueueValueType(p, e));

  } else {

    e = eqi->second;
    if ( !a->isEndPoint(p))
    {
 	m_currentPos = m_sweepLinePos;
	e->addCurveToLeft(a);
	m_currentPos = p;
	e->addCurveToRight(a);
    }
    if ( !b->isEndPoint(p) )
    {
 	m_currentPos = m_sweepLinePos;
	e->addCurveToLeft(b);
	m_currentPos = p;
	e->addCurveToRight(b);
    }
    SL_DEBUG(e->Print();)
  }
  return e;
  
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

template <class CurveInputIterator,  class SweepLineTraits_2>
inline void 
Sweep_line_base_2<CurveInputIterator, SweepLineTraits_2>::
PrintResult()
{
  std::cout << std::endl << "Result: " << std::endl;
  for ( unsigned int i = 0 ; i < m_result.size() ; i++ )
  {
    std::cout << m_result[i] << std::endl;
  }
}

CGAL_END_NAMESPACE

#endif // CGAL_SL_H
