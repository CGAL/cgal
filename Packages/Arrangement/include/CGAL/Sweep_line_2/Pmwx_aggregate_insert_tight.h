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
// release       : $CGAL_Revision: CGAL-2.3-I-81 $
// release_date  : $CGAL_Date: 2001/07/10 $
//
// file          : include/CGAL/Sweep_line_2/Pmwx_aggregate_insert_tight.h
// package       : Arrangement (2.07)
// maintainer    : Eyal Flato <flato@math.tau.ac.il>
// source        : 
// revision      : 
// revision_date : 
// author(s)     : Eti Ezra <estere@post.tau.ac.il>
//
//
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// Chapter       : 
// ======================================================================

#ifndef CGAL_PMWX_AGGREGATE_INSERT_IMPL_H
#define CGAL_PMWX_AGGREGATE_INSERT_IMPL_H

#include <list>

#include <CGAL/Sweep_line_tight_2.h>
#include <CGAL/Sweep_line_2/Pmwx_sweep_line_event.h>
#include <CGAL/Sweep_line_2/Pmwx_sweep_line_curve.h>

CGAL_BEGIN_NAMESPACE

template <class CurveInputIterator, class SweepLineTraits_2, 
          class PM_, class Change_notification_>
class Pmwx_aggregate_insert_tight :
  public Sweep_line_tight_2<CurveInputIterator, SweepLineTraits_2,
    Pmwx_sweep_line_event<SweepLineTraits_2, 
      Pmwx_sweep_line_curve<SweepLineTraits_2, 
                            typename PM_::Halfedge_handle> > ,
        Pmwx_sweep_line_curve<SweepLineTraits_2, 
                              typename PM_::Halfedge_handle> >
{
public:
  typedef SweepLineTraits_2 Traits;
  typedef typename Traits::X_curve_2 X_curve_2;
  typedef typename Traits::Point_2 Point_2;

  typedef PM_ PM;
  typedef typename PM::Halfedge_iterator  Halfedge_iterator; 
  typedef typename PM::Halfedge_handle Halfedge_handle;

  typedef Pmwx_sweep_line_curve<SweepLineTraits_2, 
                                Halfedge_handle> SubCurve;
  typedef Pmwx_sweep_line_event<SweepLineTraits_2, SubCurve> Event;
  typedef typename SubCurve::PmwxInsertInfo PmwxInsertInfo;

  typedef Change_notification_ Change_notification;
  typedef Sweep_line_tight_2<CurveInputIterator, Traits, Event, SubCurve> Base;

  typedef typename Event::VerticalXEventList VerticalXEventList;
  typedef typename Event::VerticalXEventListIter VerticalXEventListIter;

  typedef Pmwx_sweep_line_curve<SweepLineTraits_2, 
                                typename PM_::Halfedge_handle>  Subcurve;

  // repeated typedefs from the base class to avoid warnings
  typedef typename Base::EventQueueIter EventQueueIter;
  typedef typename Base::EventCurveIter EventCurveIter;
  typedef typename Base::VerticalCurveList VerticalCurveList;
  typedef typename Base::VerticalCurveListIter VerticalCurveListIter;
  typedef typename Base::StatusLineIter StatusLineIter;
  typedef typename Base::SubCurveListIter SubCurveListIter;
  typedef typename Base::SweepLinePlanarmap SweepLinePlanarmap;


  Pmwx_aggregate_insert_tight() : 
    Base(), m_change_not(NULL) {}
  
  Pmwx_aggregate_insert_tight(Traits *traits_) : 
    Base(traits_), m_change_not(NULL) {} 
  
  virtual ~Pmwx_aggregate_insert_tight() {}
  
  /*! Initializes the data structures to work with:
    - x-monotonize the input curves
    - for each end point of each curve create an event
    - for each curve in the planarmap do the same
    - initialize the event queue
    -
  */  
  void Init(CurveInputIterator begin, CurveInputIterator end, PM &pm)
  {
    Base::Init(begin, end);
    
    Halfedge_iterator eit;
    for (eit = pm.halfedges_begin(); eit != pm.halfedges_end(); ++eit, ++eit) 
    {
      InitCurve(eit->curve());
    }
    pm.clear();
  }

  void insert_curves(CurveInputIterator begin, 
                     CurveInputIterator end, 
                     PM &planarMap,
                     Change_notification* change_notification)
  {
    m_change_not = change_notification;
    std::vector<X_curve_2> subcurves;
    Init(begin, end, planarMap); 
    

    // initialize the last event in each event 
    for ( EventQueueIter qiter = m_queue->begin();
	  qiter != m_queue->end() ; ++qiter ) 
    {
      Event *e = qiter->second;
      for  (EventCurveIter rightCurveIter = e->rightCurvesBegin() ;
	    rightCurveIter != e->rightCurvesEnd() ; 
	    ++rightCurveIter )
	(*rightCurveIter)->setLastEvent(e);
      VerticalCurveList &vcurves = e->getVerticalCurves();
      VerticalCurveListIter vciter = vcurves.begin();

      while ( vciter != vcurves.end() )
      {
	if ((*vciter)->isBottomEnd(e->getPoint()))
	  (*vciter)->setLastEvent(e);
	++vciter;

      }
    }
    SL_DEBUG(
      PrintSubCurves();
      PrintEventQueue();
    )

    Sweep(planarMap, SweepLinePlanarmap());
  }
  
protected:

  /*! The main loop to calculate intersections among the curves
    Looping over the events in the queue, for each event we first
    handle the curves that are tothe left of the event point (i.e., 
    curves that we are done with), and then we look at the curves 
    to the right of the point, which means we attept to find intersections
    between them and their neighbours on the sweep line.
  */
  template <class _PM_, class Op>
  void Sweep(_PM_ &pm, Op tag)
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
      m_verticals.clear();
      m_verticalSubCurves.clear();

      while (eventIter != m_queue->end() && 
             m_traits->compare_x(eventIter->first, m_sweepLinePos) == EQUAL) {
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
	HandleLeftCurves(pm, tag);

        m_miniq.push_back(m_currentEvent);
        ++eventIter;
      }
      m_queue->erase(m_queue->begin(), eventIter);

      typename std::list<Event *>::iterator itt = m_miniq.begin();
      while ( itt != m_miniq.end())
      {
        m_currentEvent = *itt;
	HandleVerticalCurveTop(pm, tag);
        HandleRightCurves();
        ++itt;
      }
      m_miniq.clear();
      eventIter = m_queue->begin();
    }
  }

  /*! For each left-curve, if it is the "last" subcurve, i.e., the 
    event point is the right-edge of the original curve, the 
    last sub curve is created and added to the result. Otherwise
    the curve is added as is to the result.
  */
  void HandleLeftCurves(PM &pm, SweepLinePlanarmap &tag)
  {
    SL_DEBUG(std::cout << "Handling left curve" << std::endl;)
    SL_DEBUG(m_currentEvent->Print();)
    SL_DEBUG((m_currentEvent->getInsertInfo()->Print());)

    EventCurveIter leftCurveIter = m_currentEvent->leftCurvesBegin();
    m_currentPos = m_prevPos;
    const Point_2 &eventPoint = m_currentEvent->getPoint();


    Halfedge_handle h(NULL);
    m_use_hint_for_erase = false;
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
	  h = insertToPm(a, leftCurve, h, pm);
        } else {
	  h = insertToPm(cv, leftCurve, h, pm);
        }
      } else if ( leftCurve->isTarget(eventPoint))
      {
        if ( !leftCurve->isSource(lastPoint))
        {
          X_curve_2 a,b;
          m_traits->curve_split(cv, a, b, lastPoint);
	  h = insertToPm(b, leftCurve, h, pm);
        } else {
	  h = insertToPm(cv, leftCurve, h, pm);
        }

      } else { 
        X_curve_2 a,b;
        if ( leftCurve->isSource(lastPoint)) {
          m_traits->curve_split(cv, a, b, eventPoint);
	  h = insertToPm(a, leftCurve, h, pm);
        } else if ( leftCurve->isTarget(lastPoint)) {
          m_traits->curve_split(cv, b, a, eventPoint);
	  h = insertToPm(a, leftCurve, h, pm);
        } else {
          const X_curve_2 &lastCurve = leftCurve->getLastCurve();
          if ( leftCurve->isSourceLeftToTarget() ) {
            m_traits->curve_split(lastCurve, a, b, eventPoint);
	    h = insertToPm(a, leftCurve, h, pm);
          } else {
            m_traits->curve_split(lastCurve, b, a, eventPoint);
	    h = insertToPm(a, leftCurve, h, pm);
          }
        }
        leftCurve->setLastPoint(eventPoint);
        leftCurve->setLastCurve(b); 
	leftCurve->setLastEvent(m_currentEvent);
      }

      // before deleting check new neighbors that will become after deletion
      RemoveCurveFromStatusLine(leftCurve);
      m_use_hint_for_erase = true;

      m_currentPos = m_prevPos;
      ++leftCurveIter;
    }
    // when done handling the left curves, we prepare for the right curves
    m_currentEvent->initRightCurves();
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
  void HandleVerticalCurveBottom(SweepLinePlanarmap &tag)
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

      while (slIter != m_statusLine->end() &&
	     (! m_traits->curve_is_in_x_range((*slIter)->getCurve(), 
					      topEnd) ||
	      m_traits->curve_get_point_status((*slIter)->getCurve(), 
					       topEnd) != LARGER) &&
	     (! m_traits->curve_is_in_x_range((*slIter)->getCurve(), 
					      currentPoint) ||
	      m_traits->curve_get_point_status((*slIter)->getCurve(), 
					       currentPoint) != SMALLER))
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
	bool res = m_traits->nearest_intersection_to_right(cv1, cv2,
                                                           currentPoint, p, p);
	SL_DEBUG(assert(res==true);)
	  res = 0;
      
	EventQueueIter eqi = m_queue->find(p);
	Event *e = 0;
	if ( eqi == m_queue->end() )
	{
	  e = new Event(p, m_traits); 
	  SL_DEBUG(e->id = m_eventId++;)
	    m_events.push_back(e);
	  
	  e->addCurveToLeft(*slIter, m_sweepLinePos);
	  e->addCurveToRight(*slIter);
	  PRINT_NEW_EVENT(p, e);
	  m_queue->insert(EventQueueValueType(p, e));
	  
	  lastEventCreatedHere = true;

	} else {
	  e = eqi->second;
	  if ( e == prevEvent ) {
	    if ( lastEventCreatedHere )
	    {
	      if ( !(*slIter)->isLeftEnd(p) ) 
		e->addCurveToLeft(*slIter, m_sweepLinePos);
	      if ( !(*slIter)->isRightEnd(p) ) 
		e->addCurveToRight(*slIter);
	    }
	  }
	  else {
	    lastEventCreatedHere = false;
	  }
	
	  SL_DEBUG(std::cout << "Updating event \n";)
	  SL_DEBUG(e->Print();)
	}
	
	topEndEvent->addVerticalCurveXEvent(e);
	++slIter;
	prevEvent = e;
      }    
      vciter++;
    }
    
    SL_DEBUG(std::cout<<"Done Handling vertical\n";)
  }

  /*!
   *  Handle a vertical curve when the event being processed is the top end 
   *  of the curve. In this situation, the event contains a list of
   *  intersection points on the vertical curve. We go through this list and
   *  outpt the subcurves induced by these intersection points.
   *  If the curve is not vertical, returns without doing anything.
   * 
   *  @param out an iterator to the output
   *  @param tag a tag that indicates the version of the method
   */
  void HandleVerticalCurveTop(PM &pm, SweepLinePlanarmap &tag)
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

	while (slIter != m_statusLine->end() &&
	       m_traits->curve_is_in_x_range((*slIter)->getCurve(), 
					     topPoint) &&
	       m_traits->curve_get_point_status((*slIter)->getCurve(), 
						topPoint) == SMALLER &&
	       m_traits->curve_is_in_x_range((*slIter)->getCurve(), 
					     vcurve->getBottomEnd()) &&
	       m_traits->curve_get_point_status((*slIter)->getCurve(), 
					     vcurve->getBottomEnd()) == LARGER)
	{
	  SL_DEBUG(std::cout<<"checking \n";)
	  SL_DEBUG((*slIter)->Print();) 
	  if ( m_traits->compare_x((*slIter)->getLeftEnd(), topPoint) == EQUAL)
	  {
	    m_currentEvent->addVerticalCurveXEvent((*slIter)->getLastEvent(), 
						   true);
	  }
	  ++slIter;
	}   
      }

      // now we go over the list of intersection points on the vertical
      // curve in at the event and process them...
      SL_DEBUG(std::cout<<"handling the splitting now\n";)
      VerticalXEventList &pointList = m_currentEvent->getVerticalXEventList();
      if ( pointList.empty() )
      {
	insertToPmV(vcurve->getCurve(), vcurve, 
			m_currentEvent, vcurve->getLastEvent(), pm);
	++vciter;
	continue;
      }
    
      X_curve_2 a, b, c;
      a = vcurve->getCurve();
      SL_DEBUG(std::cout << "there are " << pointList.size() << " points\n";)
      Event *prevEvent = vcurve->getLastEvent();
      for ( VerticalXEventListIter i = pointList.begin() ;
	    i != pointList.end(); ++i )
      {
	SL_DEBUG(std::cout<< "splitting: " << a << " at " << *i ;)
	if ( !vcurve->isPointInRange((*i)->getPoint()) )
	{
	  SL_DEBUG(std::cout << " not !\n";)
	  continue;
	}
	SL_DEBUG(std::cout << " yes! \n";)
	m_traits->curve_split(a, b, c, (*i)->getPoint());
	if ( vcurve->isSourceLeftToTarget()) {
	  insertToPmV(b, vcurve, *i, prevEvent, pm);
	  a = c;
	} else {
	  insertToPmV(c, vcurve, *i, prevEvent, pm);
	  a = b;
	}
	vcurve->setLastEvent(*i);
	prevEvent = *i;
      }
      if ( vcurve->isSourceLeftToTarget() ) {
	insertToPmV(c, vcurve, m_currentEvent, prevEvent, pm);
      }
      else {
	insertToPmV(b, vcurve, m_currentEvent, prevEvent, pm);
      }
      ++vciter;
    }
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
  bool HandleVerticalCurveXAtEnd(Subcurve *vcurve, Subcurve *curve, 
				 Event *topEndEvent, SweepLinePlanarmap &tag)
  {
    const Point_2 &topEnd = vcurve->getTopEnd();
    // handle a curve that goes through the top point of the vertical curve
    if (m_traits->curve_is_in_x_range(curve->getCurve(), topEnd) &&
	m_traits->curve_get_point_status(curve->getCurve(), 
					 topEnd) == EQUAL)
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
    if (m_traits->curve_is_in_x_range((curve)->getCurve(), currentPoint) &&
	m_traits->curve_get_point_status((curve)->getCurve(), 
					 currentPoint) == EQUAL)
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
  
  void HandleVerticalOverlapCurves()
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
      
      if (m_traits->curve_is_in_x_range(curve->getCurve(), point) &&
	  m_traits->curve_get_point_status(curve->getCurve(), point)==SMALLER)
      {
	iter = m_verticals.erase(iter);
	
      } else if (!curve->isEndPoint(point)) {
	EventQueueIter eventIter = m_queue->find(curve->getTopEnd());
	assert(eventIter!=m_queue->end());
	(eventIter->second)->addVerticalCurveXEvent(m_currentEvent, true);
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

private:

  /*! Insert a curve to the planar map.
   *  If an identical curve was already inserted into the planarmap, it is 
   *  not inserted again. 
   *
   *  @param cv the curve to insert
   *  @param leftCurve the original curve
   *  @param hhandle a prev halfedge handle (may be NULL)
   *  @param pm a reference to the planar map
   */
  Halfedge_handle insertToPm(const X_curve_2 &cv, SubCurve *leftCurve, 
			     Halfedge_handle hhandle, PM &pm)
  {

    SL_DEBUG(std::cout << "*X inserting " << cv << "(" 
	               << leftCurve->getId() << ")\n";)
    static SubCurve *prevCurve = 0;
    static X_curve_2 prevXCv;
    
    Event *lastEvent = leftCurve->getLastEvent();
    PmwxInsertInfo *insertInfo = lastEvent->getInsertInfo();
    
    SL_DEBUG(std::cout << "lastEvent = " << lastEvent << "\n";
	     lastEvent->Print();
	     insertInfo->Print();)
    
    // if this is the same as the previous curve, don't add it again
    if ( prevCurve && SimilarCurves(cv, prevXCv)) {
      leftCurve->setLastEvent(m_currentEvent);

      if ( m_change_not ) {
	m_change_not->add_edge(cv, hhandle, true, true);
      }

      return hhandle;
    }
    prevCurve = leftCurve;
    prevXCv = cv;
    
    Halfedge_handle res; 
    
    // if the previous event on the curve is not in the planar map yet
    if ( insertInfo->getHalfedgeHandle() == Halfedge_handle(NULL) ) 
    {
      // we have a handle from the previous insert
      if ( hhandle != Halfedge_handle(NULL) ) {
	SL_DEBUG(std::cout << "  from vertex (1)";
		 std::cout << hhandle->source()->point() << " " 
		 << hhandle->target()->point() << "\n";)
        res = pm.non_intersecting_insert_from_vertex(cv, hhandle,
                                                     m_change_not);
	res = res->twin();
      } else { 
	// if this is the first left curve being inserted
	SL_DEBUG(std::cout << "  in face interior\n";)
	res = pm.insert_in_face_interior(cv, pm.unbounded_face(), m_change_not);
	if ( !leftCurve->isSourceLeftToTarget() ){
	  res = res->twin();
	}
      }
    } else 
      // the previous event on the curve is already in the planar map. 
      // Let's use it.
    {
      // skip to the right halfedge
      int jump = lastEvent->getHalfedgeJumpCount(leftCurve);
      SL_DEBUG(std::cout << "Skipping " << jump << " steps\n";)
      Halfedge_handle prev = insertInfo->getHalfedgeHandle();
      for ( int i = 0 ; i < jump ; i++ )
	prev = (prev->next_halfedge())->twin();
      
      // we have a handle from the previous insert
      if ( hhandle != Halfedge_handle(NULL) ) {
	SL_DEBUG(std::cout << "  at vertices";
		 std::cout << prev->source()->point() << " " 
		           << prev->target()->point();
		 std::cout << hhandle->source()->point() << " " 
                           << hhandle->target()->point() << "\n";)
	  res = pm.non_intersecting_insert_at_vertices(cv, prev, hhandle, 
                                                       m_change_not);
      } else {
	// if this is the first left curve being inserted
	SL_DEBUG(std::cout << "  from vertex (2)";
		 std::cout << prev->source()->point() << " " 
	                   << prev->target()->point() << "\n";)
	res = pm.non_intersecting_insert_from_vertex(cv, prev, m_change_not);
      }
    }
  
    SL_DEBUG(std::cout << "*** returning: (" << res->source()->point() << " " 
	               << res->target()->point()  << ")\n\n";)
   
    // update the information in the events so they can be used int he future
    if ( lastEvent->getNumLeftCurves() == 0 &&
	 lastEvent->isCurveLargest(leftCurve)) 
    {
      insertInfo->setHalfedgeHandle(res->twin());
    }
  
    insertInfo = m_currentEvent->getInsertInfo();
    insertInfo->setHalfedgeHandle(res);
    
    return res;

  }

  void insertToPmV(const X_curve_2 &a, SubCurve *origCurve, 
		   Event *topEvent, Event *bottomEvent, PM &pm);

  Change_notification *m_change_not;
};

/*! Insert a vertical curve to the planar map.
 *  If an identical curve was already inserted into the planarmap, it is 
 *  not inserted again. 
 *
 *  @param a the curve
 *  @param origCurve the original curve
 *  @param topEvent a pointer to the event at the top end of the curve
 *  @param bottomEvent a pointer to the event at the bottom end of the curve
 *  @param pm a reference to the planar map
 */
template <class CurveInputIterator, class SweepLineTraits_2, 
          class PM_, class Change_notification_>
void 
Pmwx_aggregate_insert_tight<CurveInputIterator, SweepLineTraits_2,
                            PM_, Change_notification_>::
insertToPmV(const X_curve_2 &cv, SubCurve *origCurve, 
	    Event *topEvent, Event *bottomEvent, PM &pm)
{
  SL_DEBUG(std::cout << "insertToPmV \n";
	   std::cout << "*V inserting " << cv << ")\n";)

  PmwxInsertInfo *topII = topEvent->getInsertInfo();
  PmwxInsertInfo *bottomII = bottomEvent->getInsertInfo();
  Halfedge_handle res; 

  // if the curve is already in the planar map, update the data and return
  if (VerticalSubCurveExists(cv)) {
    origCurve->setLastEvent(m_currentEvent);
    SL_DEBUG(std::cout << "*** returning (curve already exists\n\n";)
    return;
  }

  m_verticalSubCurves.push_back(cv);

  if ( topII->getHalfedgeHandle() == Halfedge_handle(NULL))
  {
    if ( bottomII->getHalfedgeHandle() == Halfedge_handle(NULL))
    {
      SL_DEBUG(std::cout << "  in face interior\n";)
      res = pm.insert_in_face_interior(cv, pm.unbounded_face(), m_change_not);
      if ( !origCurve->isSourceLeftToTarget() ){
	res = res->twin();
      }	
    } else 
    {
      SL_DEBUG(std::cout << "  from vertex (1) ";
	       std::cout << bottomII->getHalfedgeHandle()->source()->point() 
	                 << bottomII->getHalfedgeHandle()->target()->point()
	                 << "\n";)
      res = pm.non_intersecting_insert_from_vertex(cv, 
					  bottomII->getHalfedgeHandle(), 
						   m_change_not);
    }
  } else 
  {
    if ( bottomII->getHalfedgeHandle() == Halfedge_handle(NULL))
    {
      SL_DEBUG(std::cout << "  from vertex (2) ";
	       std::cout << topII->getHalfedgeHandle()->source()->point() 
	                 << topII->getHalfedgeHandle()->target()->point() 
                         << "\n";)
	res = pm.non_intersecting_insert_from_vertex(cv, 
					  topII->getHalfedgeHandle(),
						     m_change_not);
      res = res->twin();
    } else 
    {
      SL_DEBUG(std::cout << "  at vertices";
	       std::cout << bottomII->getHalfedgeHandle()->source()->point();
	       std::cout << bottomII->getHalfedgeHandle()->target()->point();
	       std::cout << topII->getHalfedgeHandle()->source()->point();
	       std::cout << topII->getHalfedgeHandle()->target()->point() 
                         << "\n";)
      res = pm.non_intersecting_insert_at_vertices(cv,
					     bottomII->getHalfedgeHandle(), 
					     topII->getHalfedgeHandle(),
					     m_change_not);
    }
  }

  if ( topEvent->getNumLeftCurves() == 0 ) {
    topII->setHalfedgeHandle(res);
  }

  bottomII->setHalfedgeHandle(res->twin());

  SL_DEBUG(std::cout << "*** returning: (" << res->source()->point() << " " 
	   << res->target()->point()  << ")\n\n";)

}
  



CGAL_END_NAMESPACE

#endif
