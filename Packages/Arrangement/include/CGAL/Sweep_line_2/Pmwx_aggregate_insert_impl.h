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
// file          : include/CGAL/Sweep_line_2/Pmwx_aggregate_insert_impl.h
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

#include <CGAL/Sweep_line_base_2.h>
#include <CGAL/Sweep_line_2/Pmwx_sweep_line_event.h>
#include <CGAL/Sweep_line_2/Pmwx_sweep_line_curve.h>

CGAL_BEGIN_NAMESPACE

template <class CurveInputIterator, class SweepLineTraits_2, 
          class PM_, class Change_notification_>
class Pmwx_aggregate_insert_impl :
  public Sweep_line_base_2<CurveInputIterator, SweepLineTraits_2,
    Pmwx_sweep_line_event<SweepLineTraits_2, 
      Pmwx_sweep_line_curve<SweepLineTraits_2, 
                            typename PM_::Vertex_handle, 
                            typename PM_::Halfedge_handle> > ,
        Pmwx_sweep_line_curve<SweepLineTraits_2, typename PM_::Vertex_handle, 
                              typename PM_::Halfedge_handle> >
{
public:
  typedef SweepLineTraits_2 Traits;
  typedef typename Traits::X_curve_2 X_curve_2;
  typedef typename Traits::Point_2 Point_2;

  typedef PM_ PM;
  typedef typename PM::Halfedge_iterator  Halfedge_iterator; 
  typedef typename PM::Halfedge_handle Halfedge_handle;
  typedef typename PM::Vertex_handle Vertex_handle;

  typedef Pmwx_sweep_line_curve<SweepLineTraits_2, 
                                Vertex_handle, 
                                Halfedge_handle> SubCurve;
  typedef Pmwx_sweep_line_event<SweepLineTraits_2, SubCurve> Event;
  typedef typename SubCurve::PmwxInsertInfo PmwxInsertInfo;

  typedef Change_notification_ Change_notification;
  typedef Sweep_line_base_2<CurveInputIterator, Traits, Event, SubCurve> Base;

  typedef typename Event::VerticalXEventList VerticalXEventList;
  typedef typename Event::VerticalXEventListIter VerticalXEventListIter;


  class  SweepLinePlanarmap {};

  Pmwx_aggregate_insert_impl() : 
    Base() {}
  
  Pmwx_aggregate_insert_impl(Traits *traits_) : 
    Base(traits_) {} 
  
  virtual ~Pmwx_aggregate_insert_impl() {}
  
  /*! Initializes the data structures to work with:
    - x-monotonize the inf\put curves
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
    std::vector<X_curve_2> subcurves;
    Init(begin, end, planarMap); 
    
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

  void Sweep(PM &pm, SweepLinePlanarmap tag)
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

    // intersect vertical curves...
    IntersectVerticalCurves();
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
          //*out = a; ++out;
	  h = insertToPm(a, leftCurve, h, pm);
        } else {
          //*out = cv; ++out;
	  h = insertToPm(cv, leftCurve, h, pm);
        }
      } else if ( leftCurve->isTarget(eventPoint))
      {
        if ( !leftCurve->isSource(lastPoint))
        {
          X_curve_2 a,b;
          m_traits->curve_split(cv, a, b, lastPoint);
          //*out = b; ++out;
	  h = insertToPm(b, leftCurve, h, pm);
        } else {
          //*out = cv; ++out;
	  h = insertToPm(cv, leftCurve, h, pm);
        }

      } else { 
        X_curve_2 a,b;
        if ( leftCurve->isSource(lastPoint)) {
          m_traits->curve_split(cv, a, b, eventPoint);
          //*out = a; ++out;
	  h = insertToPm(a, leftCurve, h, pm);
        } else if ( leftCurve->isTarget(lastPoint)) {
          m_traits->curve_split(cv, b, a, eventPoint);
          //*out = a; ++out;
	  h = insertToPm(a, leftCurve, h, pm);
        } else {
          const X_curve_2 &lastCurve = leftCurve->getLastCurve();
          if ( leftCurve->isSourceLeftToTarget() ) {
            m_traits->curve_split(lastCurve, a, b, eventPoint);
            //*out = a; ++out;
	    h = insertToPm(a, leftCurve, h, pm);
          } else {
            m_traits->curve_split(lastCurve, b, a, eventPoint);
            //*out = a; ++out;
	    h = insertToPm(a, leftCurve, h, pm);
          }
        }
        leftCurve->setLastPoint(eventPoint);
        leftCurve->setLastCurve(b); 
	//leftCurve->setInsertInfo(m_currentEvent->getInsertInfo());
	leftCurve->setLastEvent(m_currentEvent);
      }

      // before deleting check new neighbors that will become after deletion
      StatusLineIter sliter = 
	IntersectNeighboursAfterRemoval(leftCurve);

      m_currentPos = m_prevPos;
      m_statusLine->erase(sliter);
      ++leftCurveIter;
    }
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
	  
	} else {
	  e = eqi->second;
	  if ( !(*slIter)->isLeftEnd(p) ) 
	    e->addCurveToLeft(*slIter, m_sweepLinePos);
	  if ( !(*slIter)->isRightEnd(p) ) 
	    e->addCurveToRight(*slIter);
	  SL_DEBUG(std::cout << "Updating event \n";)
	  SL_DEBUG(e->Print();)
	}
	
	topEndEvent->addVerticalCurveXPoint(p);
	topEndEvent->addVerticalCurveXEvent(e);
	++slIter;
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
	    m_currentEvent->addVerticalCurveXEvent((*slIter)->getLastEvent(), 
						   true);
	  }
	  ++slIter;
	}   
      }

      // now we go over the list of intersection points on the vertical
      // curve in at the event and process them...
      SL_DEBUG(std::cout<<"handling the splitting now\n";)
      Halfedge_handle h(NULL);
      VerticalXEventList &pointList = m_currentEvent->getVerticalXEventList();
      if ( pointList.empty() )
      {
	h = insertToPmV(vcurve->getCurve(), vcurve, 
			m_currentEvent, vcurve->getLastEvent(), h, pm);
	++vciter;
	continue;
      }
    
      X_curve_2 a, b, c;
      a = vcurve->getCurve();
      SL_DEBUG(std::cout << "there are " << pointList.size() << " points\n";)
      SL_DEBUG(m_currentEvent->PrintVerticalXPoints();)
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
	  h = insertToPmV(b, vcurve, *i, prevEvent, h, pm);
	  a = c;
	} else {
	  h = insertToPmV(c, vcurve, *i, prevEvent, h, pm);
	  a = b;
	}
	vcurve->setLastEvent(*i);
	prevEvent = *i;
      }
      if ( vcurve->isSourceLeftToTarget() ) {
	h = insertToPmV(c, vcurve, m_currentEvent, prevEvent, h, pm);
      }
      else {
	h = insertToPmV(b, vcurve, m_currentEvent, prevEvent, h, pm);
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
      typename Traits::Curve_point_status pstatus = 
	m_traits->curve_get_point_status(curve->getCurve(), point);
      
      if ( pstatus == Traits::ABOVE_CURVE ) {
	iter = m_verticals.erase(iter);
	
      } else if (!curve->isEndPoint(point)) {
	EventQueueIter eventIter = m_queue->find(curve->getTopEnd());
	assert(eventIter!=m_queue->end());
	(eventIter->second)->addVerticalCurveXPoint(point, true);
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

  Halfedge_handle insertToPm(const X_curve_2 &a, SubCurve *leftCurve, 
			     Halfedge_handle hhandle, PM &pm)
  {
    SL_DEBUG(std::cout << "*X inserting " << a << ")\n";)

    static SubCurve *prevCurve = 0;
    static X_curve_2 prevXCv;
    Event *lastEvent = leftCurve->getLastEvent();
    SL_DEBUG(std::cout << "lastEvent = " << lastEvent << "\n";)
    assert(lastEvent!=0);
    PmwxInsertInfo *insertInfo = leftCurve->getLastEvent()->getInsertInfo();
    SL_DEBUG(insertInfo->Print();)

    if ( prevCurve && SimilarCurves(a, prevXCv)) {
      leftCurve->setLastEvent(m_currentEvent);
      return hhandle;
    }
    prevCurve = leftCurve;
    prevXCv = a;

    Halfedge_handle res; 

    if ( insertInfo->getVertexHandle() == Vertex_handle(NULL) ) {
      if ( hhandle != Halfedge_handle(NULL) ) {
	SL_DEBUG(std::cout << "  from vertex ";
		 std::cout << hhandle->target()->point() << "\n";)
	res = pm.non_intersecting_insert_from_vertex(a, hhandle, 0);
	res = res->twin();
      } else {
	SL_DEBUG(std::cout << "  in face interior\n";)
	res = pm.insert_in_face_interior(a, pm.unbounded_face(), 0);
	if ( !leftCurve->isSourceLeftToTarget() ){
	  res = res->twin();
	}
      }
    } else 
    {
      Vertex_handle v1 = insertInfo->getVertexHandle();
      if ( hhandle != Halfedge_handle(NULL) ) {
	SL_DEBUG(std::cout << "  at vertices ";
		 std::cout << v1->point();
		 std::cout << hhandle->target()->point() << "\n";)
	res = pm.non_intersecting_insert_at_vertices(a,v1,hhandle->target(),0);
      } else {
	SL_DEBUG(std::cout << "  from vertex ";
		 std::cout << v1->point() << "\n";)
	res = pm.non_intersecting_insert_from_vertex(a, v1, 0);
      }
      
    }

    SL_DEBUG(std::cout << "*** returning: (" << res->source()->point() << " " 
	     << res->target()->point()  << ")\n\n";)

    insertInfo->setVertexHandle(res->source());
    insertInfo->setHalfedgeHandle(res->twin());

    insertInfo = m_currentEvent->getInsertInfo();
    insertInfo->setVertexHandle(res->target());
    insertInfo->setHalfedgeHandle(res);

    return res;
  }
  Halfedge_handle insertToPmV(const X_curve_2 &a, SubCurve *origCurve, 
			      Event *topEvent, Event *bottomEvent,
			      Halfedge_handle hhandle, PM &pm)
  {
    SL_DEBUG(std::cout << "insertToPmV \n";)
    SL_DEBUG(std::cout << "*V inserting " << a << ")\n";)

    PmwxInsertInfo *topII = topEvent->getInsertInfo();
    PmwxInsertInfo *bottomII = bottomEvent->getInsertInfo();
    Halfedge_handle res; 

    if (VerticalSubCurveExists(a)) {
      origCurve->setLastEvent(m_currentEvent);
      SL_DEBUG(std::cout << "*** returning (curve already exists\n\n";)
      return hhandle;
    }
    m_verticalSubCurves.push_back(a);

    if ( topII->getVertexHandle() == Vertex_handle(NULL))
    {
      if ( bottomII->getVertexHandle() == Vertex_handle(NULL))
      {
	SL_DEBUG(std::cout << "  in face interior\n";)
	res = pm.insert_in_face_interior(a, pm.unbounded_face(), 0);
	if ( !origCurve->isSourceLeftToTarget() ){
	  res = res->twin();
	}	
      } else 
      {
	SL_DEBUG(std::cout << "  from vertex (2)";
		 std::cout << bottomII->getVertexHandle()->point() << "\n";)
	res = pm.non_intersecting_insert_from_vertex(a, 
					     bottomII->getVertexHandle(), 0);
      }
    } else 
    {
      if ( bottomII->getVertexHandle() == Vertex_handle(NULL))
      {
	SL_DEBUG(std::cout << "  from vertex (2)";
		 std::cout << topII->getVertexHandle()->point() << "\n";)
	res = pm.non_intersecting_insert_from_vertex(a, 
					     topII->getVertexHandle(), 0);
	res = res->twin();
      } else 
      {
	SL_DEBUG(std::cout << "  at vertices ";
		 std::cout << bottomII->getVertexHandle()->point();
		 std::cout << topII->getVertexHandle()->point() << "\n";)
	res = pm.non_intersecting_insert_at_vertices(a,
					     bottomII->getVertexHandle(), 
					     topII->getVertexHandle(),0);
      }
    }

    topII->setVertexHandle(res->target());
    topII->setHalfedgeHandle(res);

    bottomII->setVertexHandle(res->source());
    bottomII->setHalfedgeHandle(res->twin());

    SL_DEBUG(std::cout << "*** returning: (" << res->source()->point() << " " 
	     << res->target()->point()  << ")\n\n";)

    return res;
  }
  
};

CGAL_END_NAMESPACE

#endif
