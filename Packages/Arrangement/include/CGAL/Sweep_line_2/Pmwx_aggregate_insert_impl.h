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

#define A(x)

CGAL_BEGIN_NAMESPACE

template <class CurveInputIterator, class SweepLineTraits_2, 
          class PM_, class Change_notification_>
class Pmwx_aggregate_insert_impl :
  public Sweep_line_base_2<CurveInputIterator, SweepLineTraits_2,
       Pmwx_sweep_line_event<SweepLineTraits_2, 
                             Pmwx_sweep_line_curve<SweepLineTraits_2, 
                                              typename PM_::Vertex_handle, 
                                              typename PM_::Halfedge_handle> > ,
       Pmwx_sweep_line_curve<SweepLineTraits_2, typename  PM_::Vertex_handle, 
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
	//HandleVerticalCurveBottom(tag);
        HandleLeftCurves(pm, tag);

        m_miniq.push_back(m_currentEvent);
        ++eventIter;
      }
      m_queue->erase(m_queue->begin(), eventIter);

      typename std::list<Event *>::iterator itt = m_miniq.begin();
      while ( itt != m_miniq.end())
      {
        m_currentEvent = *itt;
	//HandleVerticalCurveTop(out, tag);
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

      A(std::cout << "\n\n---- " << m_currentEvent->getPoint() << "----- \n";)
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
	leftCurve->setInsertInfo(m_currentEvent->getInsertInfo());
      }

      // before deleting check new neighbors that will become after deletion
      StatusLineIter sliter = 
	IntersectNeighboursAfterRemoval(leftCurve);

      m_currentPos = m_prevPos;
      m_statusLine->erase(sliter);
      ++leftCurveIter;
    }
  }

private:

  Halfedge_handle insertToPm(const X_curve_2 &a, SubCurve *leftCurve, 
			     Halfedge_handle hhandle, PM &pm)
  {

    A(std::cout << "inserting " << a << "(" << leftCurve->getId() << ")\n";)
    Halfedge_handle res; 
    PmwxInsertInfo *insertInfo = leftCurve->getInsertInfo();

    if ( insertInfo == 0 ) {
      A(std::cout << "v1 is out\n";)
      if ( hhandle != Halfedge_handle(NULL) ) {
	A(std::cout << "  from vertex\n";
	  std::cout << hhandle->target()->point() << "\n";)
	res = pm.non_intersecting_insert_from_vertex(a, hhandle, 0);
	res = res->twin();
      } else {
	A(std::cout << "  insert\n";)
	  //res = pm.non_intersecting_insert(a, 0);
	res = pm.insert_in_face_interior(a, pm.unbounded_face(), 0);
	if ( !leftCurve->isSourceLeftToTarget() ){
	  res = res->twin();
	}

      }
    } else 
    {
      A(std::cout << "v1 is in\n";)
      Vertex_handle v1 = insertInfo->getVertexHandle();
      A(std::cout << "(" << v1->point() <<  ")\n";)
      if ( hhandle != Halfedge_handle(NULL) ) {
	//res = pm.non_intersecting_insert_from_vertex(a, hhandle, 0);
	A(std::cout << "  at vertices\n";
	  std::cout << v1->point() << "\n";
	  std::cout << hhandle->target()->point() << "\n";)
	res = pm.non_intersecting_insert_at_vertices(a,v1,hhandle->target(),0);
      } else {
	A(std::cout << "  from vertex\n";
	  std::cout << v1->point() << "\n";)
	res = pm.non_intersecting_insert_from_vertex(a, v1, 0);
	//res = res->twin();
      }
      
    }

    A(std::cout << "\n*** returning: (" << res->source()->point() << " " 
      << res->target()->point()  << ")\n";)

    insertInfo = m_currentEvent->getInsertInfo();
    insertInfo->setVertexHandle(res->target());
    insertInfo->setHalfedgeHandle(res);

    return res;
  }

};
CGAL_END_NAMESPACE

#endif
