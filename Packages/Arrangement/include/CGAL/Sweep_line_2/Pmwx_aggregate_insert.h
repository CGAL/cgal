// Copyright (c) 1997  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Tali Zvi <talizvi@post.tau.ac.il>

#ifndef CGAL_PMWX_AGGREGATE_INSERT_H
#define CGAL_PMWX_AGGREGATE_INSERT_H

#include <CGAL/Sweep_line_2/Sweep_line_2_impl.h>
#include <CGAL/Sweep_line_2/Pmwx_sweep_line_event.h>
#include <CGAL/Sweep_line_2/Pmwx_sweep_line_curve.h>
#include <CGAL/assertions.h>
#include <list>

CGAL_BEGIN_NAMESPACE

template <class CurveInputIterator, class SweepLineTraits_2, 
          class PM_, class Change_notification_>
class Pmwx_aggregate_insert :
  public Sweep_line_2_impl<CurveInputIterator, SweepLineTraits_2,
    Pmwx_sweep_line_event<SweepLineTraits_2, 
      Pmwx_sweep_line_curve<SweepLineTraits_2, 
                            typename PM_::Halfedge_handle> > ,
        Pmwx_sweep_line_curve<SweepLineTraits_2, 
                              typename PM_::Halfedge_handle> >
{
public:
  typedef SweepLineTraits_2 Traits;
  typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;
  typedef typename Traits::Point_2 Point_2;

  typedef PM_ PM;
  typedef typename PM::Halfedge_iterator  Halfedge_iterator; 
  typedef typename PM::Halfedge_handle Halfedge_handle;

  

  typedef Pmwx_sweep_line_curve<SweepLineTraits_2, 
                                Halfedge_handle> SubCurve;
  typedef Pmwx_sweep_line_event<SweepLineTraits_2, SubCurve> Event;
  typedef typename SubCurve::PmwxInsertInfo PmwxInsertInfo;

  typedef Change_notification_ Change_notification;
  typedef Sweep_line_2_impl<CurveInputIterator, Traits, Event, SubCurve> Base;

  typedef typename Event::VerticalXEventSet VerticalXEventSet;
  typedef typename Event::VerticalXEventSetIter VerticalXEventSetIter;

  typedef Pmwx_sweep_line_curve<SweepLineTraits_2, 
                                typename PM_::Halfedge_handle>  Subcurve;

  // repeated typedefs from the base class to avoid warnings
  typedef typename Base::EventQueueIter EventQueueIter;
  typedef typename Base::EventCurveIter EventCurveIter;
  //typedef typename Base::VerticalCurveList VerticalCurveList;
  //typedef typename Base::VerticalCurveListIter VerticalCurveListIter;
  typedef typename Base::StatusLineIter StatusLineIter;
  typedef typename Base::SubCurveList SubCurveList;
  typedef typename Base::SubCurveListIter SubCurveListIter;
  typedef typename Base::SweepLinePlanarmap SweepLinePlanarmap;
  typedef typename Base::EventQueueValueType EventQueueValueType;

#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_3
  using Base::m_queue;
  //using Base::m_prevPos;
  using Base::m_traits;
  using Base::m_sweepLinePos;
  //using Base::m_verticals;
  //using Base::m_verticalSubCurves;
 // using Base::m_currentPos;
  using Base::m_currentEvent;
  //using Base::m_miniq;
  using Base::m_use_hint_for_erase;
  using Base::m_statusLine;
  using Base::m_tmpOut;
  using Base::m_status_line_insert_hint;
  //using Base::m_events;
#endif

  Pmwx_aggregate_insert() : 
    Base(), m_change_not(NULL)
    {
    }
  
  Pmwx_aggregate_insert(Traits *traits_) : 
    Base(traits_), m_change_not(NULL) {} 
  
  virtual ~Pmwx_aggregate_insert() {}
  
  /*! Initializes the data structures to work with:
    - x-monotonize the input curves
    - for each end point of each curve create an event
    - for each curve in the planarmap do the same
    - initialize the event queue
    -
  */  
 /* void init(CurveInputIterator begin, CurveInputIterator end, PM &pm)
  {
    Base::init(begin, end);
    Halfedge_iterator eit;

    for (eit = pm.halfedges_begin(); eit != pm.halfedges_end(); ++eit, ++eit) 
    {
      init_curve(eit->curve());
    }
    pm.clear();
  }*/

  void init(CurveInputIterator begin, CurveInputIterator end, PM &pm)
  {
   // Base::init(begin, end);
    Halfedge_iterator eit;

    for (eit = pm.halfedges_begin(); eit != pm.halfedges_end(); ++eit, ++eit) 
    {
      //init_curve(eit->curve());
      m_xcurves.push_back(eit->curve());
      ++m_num_of_subCurves;
    }
    pm.clear();
    Base::init(begin, end);
  }

  void insert_curves(CurveInputIterator begin, 
                     CurveInputIterator end, 
                     PM &planarMap,
                     Change_notification* change_notification)
  {
    m_change_not = change_notification;
    std::vector<X_monotone_curve_2> subcurves;
    init(begin, end, planarMap); 
    

    // initialize the last event in each event 
    for ( EventQueueIter qiter = m_queue->begin();
          qiter != m_queue->end() ; ++qiter ) 
    {
      Event *e = qiter->second;
      for  (EventCurveIter rightCurveIter = e->right_curves_begin() ;
            rightCurveIter != e->right_curves_end() ; 
            ++rightCurveIter )
        (*rightCurveIter)->set_last_event(e);
     /* VerticalCurveList &vcurves = e->get_vertical_curves();
      VerticalCurveListIter vciter = vcurves.begin();

      while ( vciter != vcurves.end() )
      {
        if ((*vciter)->is_bottom_end(e->get_point()))
          (*vciter)->set_last_event(e);
        ++vciter;

      }*/
    }
    SL_DEBUG(
      PrintSubCurves();
      PrintEventQueue();
      );

    sweep(planarMap, SweepLinePlanarmap());
  }
  
protected:

  /*! The main loop to calculate intersections among the curves
    Looping over the events in the queue, for each event we first
    handle the curves that are to the left of the event point (i.e., 
    curves that we are done with), and then we look at the curves 
    to the right of the point, which means we attempt to find intersections
    between them and their neighbours on the sweep line.
  */
  template <class _PM_, class Op>
  void sweep(_PM_ &pm, Op tag)
  {
    EventQueueIter eventIter = m_queue->begin();
    //m_prevPos = eventIter->first;
   // Point_2 referencePoint;

    while ( eventIter != m_queue->end() )
    {
      const Point_2 &p = (*eventIter).first;
      //if ( m_traits->compare_x(m_sweepLinePos, p) == SMALLER ) 
      //{
      // // m_prevPos = m_sweepLinePos;
      //  m_verticals.clear();
      //  m_verticalSubCurves.clear();
      //}
      m_sweepLinePos = p;
     // m_currentPos = p;

     // p = (*eventIter).first;
      m_currentEvent = eventIter->second;
      SL_DEBUG(std::cout << "------------- " << p << " --------------"
                         << std::endl;
               PrintStatusLine();
               m_currentEvent->Print();
               );

      /*if ( m_traits->compare_x((*eventIter).first, m_sweepLinePos) !=
           EQUAL)
      {
        SL_DEBUG(std::cout << "clearing miniq " 
                           << (*eventIter).first  << " "
                           << "\n";);
        m_miniq.clear();
      }
      m_miniq.push_back(m_currentEvent);*/

      //handle_vertical_curve_bottom(tag);
      //handle_vertical_overlap_curves();
      handle_left_curves(pm, tag);

      m_queue->erase(eventIter);
     
      //handle_vertical_curve_top(pm, tag);
      handle_right_curves(pm, tag);
      if(m_currentEvent->get_insert_info()->get_right_curves_counter() == 0)
      {
        //delete m_currentEvent;
        //m_events.remove(m_currentEvent);
         m_eventAlloc.destroy(m_currentEvent);
         m_eventAlloc.deallocate(m_currentEvent,1);
      }
        
      
      eventIter = m_queue->begin();
    }
  }

  /*! For each left-curve, if it is the "last" subcurve, i.e., the 
    event point is the right-edge of the original curve, the 
    last sub curve is created and added to the result. Otherwise
    the curve is added as is to the result.
  */
  void handle_left_curves(PM &pm, SweepLinePlanarmap &tag)
  {
    SL_DEBUG(std::cout << "Handling left curve" << std::endl;);
    SL_DEBUG(m_currentEvent->Print(););
    SL_DEBUG((m_currentEvent->get_insert_info()->Print()););

   if(m_currentEvent->does_contain_overlap())
     handle_left_curves_overlap(pm, tag);
   else
     handle_left_curves_no_overlap(pm, tag);
  }



  ///////////////////////////////////////////////////////////////////////////
  void handle_left_curves_overlap(PM &pm, SweepLinePlanarmap &tag)
  {

    EventCurveIter leftCurveIter = m_currentEvent->left_curves_begin();
    //m_currentPos = m_prevPos;
    const Point_2 &eventPoint = m_currentEvent->get_point();

    Halfedge_handle h(NULL);
    m_use_hint_for_erase = false;
    bool are_overlap = false;
    SubCurve* leftCurvePrev = 0;
    while ( leftCurveIter != m_currentEvent->left_curves_end() )  
    {
      SubCurve *leftCurve = *leftCurveIter;
      if(leftCurvePrev && do_curves_overlap(leftCurvePrev,leftCurve) )
         are_overlap = true;

       
      const X_monotone_curve_2 &cv = leftCurve->get_curve();
      const Point_2 &lastPoint = leftCurve->get_last_point();
      
      if ( leftCurve->is_source(eventPoint))
      {
        if ( !leftCurve->is_target(lastPoint) )
        { 
          X_monotone_curve_2 a ,b;
          m_traits->curve_split(cv, a, b, lastPoint);
          if(!are_overlap)
            h = insert_to_pm(a, leftCurve, h, pm);
          else
             if(leftCurve->get_last_event()->get_insert_info()->dec_right_curves_counter() == 0)
             {
               //delete leftCurve->get_last_event();
                 m_eventAlloc.destroy(leftCurve->get_last_event());
                 m_eventAlloc.deallocate(leftCurve->get_last_event(),1);
             }
        }
        else 
        {
          if(!are_overlap)
            h = insert_to_pm(cv, leftCurve, h, pm);
          else
             if(leftCurve->get_last_event()->get_insert_info()->dec_right_curves_counter() == 0)
             {
               //delete leftCurve->get_last_event();
                 m_eventAlloc.destroy(leftCurve->get_last_event());
                 m_eventAlloc.deallocate(leftCurve->get_last_event(),1);
             }
        }
      }
      else 
        if ( leftCurve->is_target(eventPoint))
        {
          if ( !leftCurve->is_source(lastPoint))
          {
            X_monotone_curve_2 a ,b;
            m_traits->curve_split(cv, a, b, lastPoint);
            if(!are_overlap)
              h = insert_to_pm(b, leftCurve, h, pm);
            else
             if(leftCurve->get_last_event()->get_insert_info()->dec_right_curves_counter() == 0)
             {
               //delete leftCurve->get_last_event();
                 m_eventAlloc.destroy(leftCurve->get_last_event());
                 m_eventAlloc.deallocate(leftCurve->get_last_event(),1);
             }
          }
          else 
          {
            if(!are_overlap)
              h = insert_to_pm(cv, leftCurve, h, pm);
            else
             if(leftCurve->get_last_event()->get_insert_info()->dec_right_curves_counter() == 0)
             {
               //delete leftCurve->get_last_event();
                 m_eventAlloc.destroy(leftCurve->get_last_event());
                 m_eventAlloc.deallocate(leftCurve->get_last_event(),1);
             }
          }
        }
        else  // the event point passes through the interior of 'cv'
        { 
          X_monotone_curve_2 a ,b;
          if ( leftCurve->is_source(lastPoint))
          {
            m_traits->curve_split(cv, a, b, eventPoint);
            if(!are_overlap)
              h = insert_to_pm(a, leftCurve, h, pm);
            else
             if(leftCurve->get_last_event()->get_insert_info()->dec_right_curves_counter() == 0)
             {
               //delete leftCurve->get_last_event();
                 m_eventAlloc.destroy(leftCurve->get_last_event());
                 m_eventAlloc.deallocate(leftCurve->get_last_event(),1);
             }
          }  
          else
            if ( leftCurve->is_target(lastPoint))
            {
              m_traits->curve_split(cv, b, a, eventPoint);
              if(!are_overlap)
                h = insert_to_pm(a, leftCurve, h, pm);
              else
               if(leftCurve->get_last_event()->get_insert_info()->dec_right_curves_counter() == 0)
               {
                 //delete leftCurve->get_last_event();
                   m_eventAlloc.destroy(leftCurve->get_last_event());
                   m_eventAlloc.deallocate(leftCurve->get_last_event(),1);
               }
            }
            else 
            {
              const X_monotone_curve_2 &lastCurve = leftCurve->get_last_curve();
              if ( leftCurve->is_source_left_to_target() )
              {
                m_traits->curve_split(lastCurve, a, b, eventPoint);
                if(!are_overlap)
                  h = insert_to_pm(a, leftCurve, h, pm);
                else
                 if(leftCurve->get_last_event()->get_insert_info()->dec_right_curves_counter() == 0)
                 {
                   //delete leftCurve->get_last_event();
                    m_eventAlloc.destroy(leftCurve->get_last_event());
                    m_eventAlloc.deallocate(leftCurve->get_last_event(),1);
                 }
              }
              else 
              {
                m_traits->curve_split(lastCurve, b, a, eventPoint);
                if(!are_overlap)
                  h = insert_to_pm(a, leftCurve, h, pm);
                else
                 if(leftCurve->get_last_event()->get_insert_info()->dec_right_curves_counter() == 0)
                 {
                   //delete leftCurve->get_last_event();
                    m_eventAlloc.destroy(leftCurve->get_last_event());
                    m_eventAlloc.deallocate(leftCurve->get_last_event(),1);
                 } 
              }
            }
            leftCurve->set_last_point(eventPoint);
            leftCurve->set_last_curve(b); 
            leftCurve->set_last_event(m_currentEvent);
        }

        // before deleting check new neighbors that will become after deletion
        remove_curve_from_status_line(leftCurve);
        m_use_hint_for_erase = true;
      
      //  m_currentPos = m_prevPos;
        leftCurvePrev = *leftCurveIter; 
        ++leftCurveIter;
        are_overlap = false;
        
    }  
    
    // when done handling the left curves, we prepare for the right curves
    m_currentEvent->init_right_curves();
  }
  /////////////////////////////////////////////////////////////////////////////
  void handle_left_curves_no_overlap(PM &pm, SweepLinePlanarmap &tag)
  {
    EventCurveIter leftCurveIter = m_currentEvent->left_curves_begin();
    //m_currentPos = m_prevPos;
    const Point_2 &eventPoint = m_currentEvent->get_point();

    Halfedge_handle h(NULL);
    m_use_hint_for_erase = false;
    while ( leftCurveIter != m_currentEvent->left_curves_end() )  
    {
      SubCurve *leftCurve = *leftCurveIter;   
      const X_monotone_curve_2 &cv = leftCurve->get_curve();
      const Point_2 &lastPoint = leftCurve->get_last_point();
      
      if ( leftCurve->is_source(eventPoint))
      {
        if ( !leftCurve->is_target(lastPoint) )
        { 
          X_monotone_curve_2 a ,b;
          m_traits->curve_split(cv, a, b, lastPoint);
          h = insert_to_pm(a, leftCurve, h, pm);
        }
        else 
        {
          h = insert_to_pm(cv, leftCurve, h, pm);
        }
      }
      else 
        if ( leftCurve->is_target(eventPoint))
        {
          if ( !leftCurve->is_source(lastPoint))
          {
            X_monotone_curve_2 a ,b;
            m_traits->curve_split(cv, a, b, lastPoint);
            h = insert_to_pm(b, leftCurve, h, pm);
          }
          else 
          {
            h = insert_to_pm(cv, leftCurve, h, pm);
          }
        }
        else  // the event point passes through the interior of 'cv'
        { 
          X_monotone_curve_2 a ,b;
          if ( leftCurve->is_source(lastPoint))
          {
            m_traits->curve_split(cv, a, b, eventPoint);
            h = insert_to_pm(a, leftCurve, h, pm);
          }  
          else
            if ( leftCurve->is_target(lastPoint))
            {
              m_traits->curve_split(cv, b, a, eventPoint);
              h = insert_to_pm(a, leftCurve, h, pm);
            }
            else 
            {
              const X_monotone_curve_2 &lastCurve = leftCurve->get_last_curve();
              if ( leftCurve->is_source_left_to_target() )
              {
                m_traits->curve_split(lastCurve, a, b, eventPoint);
                h = insert_to_pm(a, leftCurve, h, pm);
              }
              else 
              {
                m_traits->curve_split(lastCurve, b, a, eventPoint);
                h = insert_to_pm(a, leftCurve, h, pm);
                
              }
            }
            leftCurve->set_last_point(eventPoint);
            leftCurve->set_last_curve(b); 
            leftCurve->set_last_event(m_currentEvent);
        }

        // before deleting check new neighbors that will become after deletion
        remove_curve_from_status_line(leftCurve);
        m_use_hint_for_erase = true;
      
      //  m_currentPos = m_prevPos;
        ++leftCurveIter;  
    }  
    
    // when done handling the left curves, we prepare for the right curves
    m_currentEvent->init_right_curves();
  }
  ////////////////////////////////////////////////////////////////////////


 
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
  //void handle_vertical_curve_bottom(SweepLinePlanarmap &tag)
  //{
  //  SL_DEBUG(std::cout << "\nhandle_vertical_curve_bottom... ("
  //                     << m_currentEvent->get_point() << ")\n";);
  //  if ( !m_currentEvent->does_contain_vertical_curve() )
  //  {
  //    SL_DEBUG(std::cout<<" - not vertical - exiting\n ";);
  //    return;
  //  } 
  //  SL_DEBUG(std::cout<<"\n ";);

  //  VerticalCurveList &vcurves = m_currentEvent->get_vertical_curves();
  //  VerticalCurveListIter vciter = vcurves.begin();
  //  const Point_2 &currentPoint = m_currentEvent->get_point();
  //  
  //  SL_DEBUG(std::cout << vcurves.size() << " vertical curves in event\n";);
  //  while ( vciter != vcurves.end() )
  //  {
  //    Subcurve *vcurve = *vciter;
  //    SL_DEBUG(std::cout << "working on " << vcurve->get_curve() << "\n";);
  //    if ( vcurve->is_top_end(currentPoint))
  //    {
  //      vciter++;
  //      continue;
  //    }
  //  
  //    SL_DEBUG(std::cout<<"handling bottom point of vertical curve\n";);
  //    StatusLineIter slIter = m_statusLine->lower_bound(vcurve);
  //    if ( slIter == m_statusLine->end() ) 
  //    {
  //      SL_DEBUG(std::cout<<"no curves intersecting. exiting\n";);
  //      vciter++;
  //      continue;
  //    }    
  //    
  //    SL_DEBUG(std::cout<<"starting at curve \n";);
  //    SL_DEBUG((*slIter)->Print(););
  //    const Point_2 &topEnd = vcurve->get_top_end();
  //    EventQueueIter topEndEventIter = m_queue->find(topEnd);
  //    CGAL_assertion(topEndEventIter != m_queue->end());
  //    Event *topEndEvent = topEndEventIter->second;
  //    
  //    bool lastEventCreatedHere = false;
  //    Event *prevEvent = 0;

  //    while (slIter != m_statusLine->end() &&
  //           (! m_traits->point_in_x_range((*slIter)->get_curve(), 
  //                                            topEnd) ||
  //            m_traits->curve_compare_y_at_x(topEnd, (*slIter)->get_curve()) !=
  //            SMALLER) &&
  //           (! m_traits->point_in_x_range((*slIter)->get_curve(), 
  //                                            currentPoint) ||
  //            m_traits->curve_compare_y_at_x(currentPoint,
  //                                       (*slIter)->get_curve()) != LARGER))
  //    {
  //      SL_DEBUG(std::cout<<"intersecting with \n";);
  //      SL_DEBUG((*slIter)->Print(););
  //        
  //      if ( handle_vertical_curve_x_at_end(vcurve, *slIter, topEndEvent, tag))
  //      {
  //        ++slIter;
  //        continue;
  //      }
  //      
  //      // handle a curve that goes through the interior of the vertical curve
  //      const X_monotone_curve_2 &cv1 = vcurve->get_curve();
  //      const X_monotone_curve_2 &cv2 = (*slIter)->get_curve();
  //      Object res = m_traits->nearest_intersection_to_right(cv1, cv2,
  //                                                           currentPoint);
  //      Point_2 p;
  //      //Point_2 *p = new Point_2(); 
  //      // m_points.push_back(p);

  //      if (!CGAL::assign(p, res))
  //        CGAL_assertion(0);
  //    
  //       const std::pair<EventQueueIter,bool>& insert_res = 
  //        (m_queue->insert(EventQueueValueType(p,0)));

  //      Event *e ;
  //      if(insert_res.second)    
  //      {
  //        e = new Event(p, m_traits); 
  //        SL_DEBUG(e->id = m_eventId++;);
  //         // m_events.push_back(e); 
  //        
  //        e->add_curve_to_left(*slIter/*, m_sweepLinePos*/);
  //        e->add_curve_to_right(*slIter);
  //      //  e->mark_internal_intersection_point(); //Baruch
  //        PRINT_NEW_EVENT(p, e);
  //        //m_queue->insert(EventQueueValueType(&(e->get_point()), e));

  //        
  //        lastEventCreatedHere = true;
  //        (insert_res.first)->second = e;


  //      } else {
  //       // e = eqi->second;
  //        e = (insert_res.first)->second;

  //        if ( e == prevEvent ) {
  //          if ( lastEventCreatedHere )
  //          {
  //            bool isLeftEnd = (*slIter)->is_left_end(p) ,
  //                 isRightEnd = (*slIter)->is_right_end(p);
  //            if( !isLeftEnd && !isRightEnd)
  //            {
  //              e->add_curve_to_left(*slIter/*, m_sweepLinePos*/);
  //              e->add_curve_to_right(*slIter);
  //             // e->mark_internal_intersection_point(); //Baruch
  //            }
  //            else
  //              if ( !(*slIter)->is_left_end(p) ) 
  //                e->add_curve_to_left(*slIter/*, m_sweepLinePos*/);
  //              else
  //                if ( !(*slIter)->is_right_end(p) ) 
  //                  e->add_curve_to_right(*slIter);
  //           }
  //        }
  //        else {
  //          lastEventCreatedHere = false;
  //        }
  //      
  //        SL_DEBUG(std::cout << "Updating event \n";);
  //        SL_DEBUG(e->Print(););
  //      }
  //      
  //      topEndEvent->add_vertical_curve_x_event(e);
  //      ++slIter;
  //      prevEvent = e;
  //    }    
  //    vciter++;
  //  }
  //  
  //  SL_DEBUG(std::cout<<"Done Handling vertical\n";);
  //}

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
  //void handle_vertical_curve_top(PM &pm, SweepLinePlanarmap &tag)
  //{
  //  SL_DEBUG(std::cout << "handle_vertical_curve_top... (" 
  //                     << m_currentEvent->get_point() << ")\n";);
  //  if ( !m_currentEvent->does_contain_vertical_curve() ) {
  //    SL_DEBUG(std::cout<<"exiting\n ";);
  //    return;
  //  }
  //  SL_DEBUG(std::cout<<"\n ";);

  //  VerticalCurveList &vcurves = m_currentEvent->get_vertical_curves();
  //  VerticalCurveListIter vciter = vcurves.begin();

  //  while ( vciter !=vcurves.end() )
  //  {

  //    Subcurve *vcurve = *vciter;
  //    const Point_2 &topPoint = m_currentEvent->get_point();
  //    // if this is the bottom point, nothing to do here
  //    if ( vcurve->is_bottom_end(topPoint)) {
  //      SL_DEBUG(std::cout<<"this is the bottom. skipping.\n";);
  //      ++vciter;
  //      continue;
  //    }

  //    SL_DEBUG(std::cout<<"handling top point of vertical curve\n";);


  //    // the following while loop comes to handle  | 
  //    // in the case where a new curve begins at   |------
  //    // a vertical curve                          |

  //    // find the "position" of the curve of the status line
  //    StatusLineIter slIter = m_statusLine->lower_bound(vcurve);
  //    
  //    if ( slIter != m_statusLine->end() ) 
  //    {
  //      SL_DEBUG(std::cout<<"starting at curve \n";);
  //      SL_DEBUG((*slIter)->Print(););



  //      for( ; slIter != m_statusLine->end() ; ++slIter)
  //      {
  //        if(  m_traits->point_in_x_range((*slIter)->get_curve(), 
  //                                           topPoint) &&
  //             m_traits->curve_compare_y_at_x(topPoint,
  //                                            (*slIter)->get_curve()) ==
  //             LARGER &&
  //             m_traits->point_in_x_range((*slIter)->get_curve(), 
  //                                           vcurve->get_bottom_end()) &&
  //             m_traits->curve_compare_y_at_x(vcurve->get_bottom_end(),
  //                                            (*slIter)->get_curve()) == 
  //             SMALLER)
  //        {
  //          SL_DEBUG(std::cout<<"checking \n";);
  //          SL_DEBUG((*slIter)->Print(););
  //          if ( m_traits->compare_x((*slIter)->get_left_end(), topPoint) == 
  //                                                                     EQUAL)
  //          {
  //            m_currentEvent->add_vertical_curve_x_event(
  //                                             (*slIter)->get_last_event());
  //                                             
  //          }
  //        }
  //      }

  //      //  (comment the next whilw loop)
  //      /*while (slIter != m_statusLine->end() &&
  //             m_traits->point_in_x_range((*slIter)->get_curve(), 
  //                                           topPoint) &&
  //             m_traits->curve_compare_y_at_x(topPoint,
  //                                            (*slIter)->get_curve()) ==
  //             LARGER &&
  //             m_traits->point_in_x_range((*slIter)->get_curve(), 
  //                                           vcurve->get_bottom_end()) &&
  //             m_traits->curve_compare_y_at_x(vcurve->get_bottom_end(),
  //                                            (*slIter)->get_curve()) == 
  //             SMALLER)
  //      {
  //        SL_DEBUG(std::cout<<"checking \n";);
  //        SL_DEBUG((*slIter)->Print(););
  //        if ( m_traits->compare_x((*slIter)->get_left_end(), topPoint) == 
  //                                                                     EQUAL)
  //        {
  //          m_currentEvent->add_vertical_curve_x_event(
  //                                             (*slIter)->get_last_event(), 
  //                                             true);
  //        }
  //        ++slIter;
  //      }   */
  //    }

  //    // now we go over the list of intersection points on the vertical
  //    // curve in at the event and process them...
  //    SL_DEBUG(std::cout<<"handling the splitting now\n";);
  //    VerticalXEventSet &pointList = 
  //                           m_currentEvent->get_vertical_x_event_set();
  //    if ( pointList.empty() )
  //    {
  //      //in that case there aren't any intersection point on the vertical
  //      //curve , so no need to split at. just insert it to the planar map
  //      insert_to_pm_v(vcurve->get_curve(), vcurve, 
  //                      m_currentEvent, vcurve->get_last_event(), pm);
  //      ++vciter;
  //      continue;
  //    }
  //  
  //    X_monotone_curve_2 a (vcurve->get_curve());
  //    SL_DEBUG(std::cout << "there are " << pointList.size() << " points\n";);
  //    Event *prevEvent = vcurve->get_last_event();
  //    for ( VerticalXEventSetIter i = pointList.begin() ;
  //          i != pointList.end(); ++i )
  //    {
  //      SL_DEBUG(std::cout<< "splitting: " << a << " at " <<(*i)->get_point(););
  //      if ( !vcurve->is_point_in_range((*i)->get_point()) )
  //      {
  //        SL_DEBUG(std::cout << " not !\n";);
  //        continue;
  //      }
  //      SL_DEBUG(std::cout << " yes! \n";);

  //      X_monotone_curve_2  b, c;
  //      m_traits->curve_split(a, b, c, (*i)->get_point());
  //      if ( vcurve->is_source_left_to_target()) {
  //        insert_to_pm_v(b, vcurve, *i, prevEvent, pm);
  //        a = c;
  //      } else {
  //        insert_to_pm_v(c, vcurve, *i, prevEvent, pm);
  //        a = b;
  //      }
  //      vcurve->set_last_event(*i);
  //      prevEvent = *i;
  //    }
  //  
  // 
  //    insert_to_pm_v(a,vcurve, m_currentEvent , prevEvent, pm); 
  // 
  //    ++vciter;
  //  }
  //}

  /*!
   * Handles the case in which a curve on the status line passes through
   * one of the end points of the vertical curve.
   *
   * @param vcurve the vertical curve we are dealing with
   * @param curve a cerve that intersects with the vertical curve
   * @param topEndEvent the event attached to the top end of the vertical curve
   * @param tag 
   * @return returns true if the curve passed through one of the ends of the 
   *              vertical curve. Returns false otherwise.
   */
  //bool handle_vertical_curve_x_at_end(Subcurve *vcurve, Subcurve *curve, 
  //                               Event *topEndEvent, SweepLinePlanarmap &tag)
  //{
  //  const Point_2 &topEnd = vcurve->get_top_end();
  //  // handle a curve that goes through the top point of the vertical curve
  //  if (m_traits->point_in_x_range(curve->get_curve(), topEnd) &&
  //      m_traits->curve_compare_y_at_x(topEnd, curve->get_curve()) == EQUAL)
  //  {
  //    bool isLeftEnd = curve->is_left_end(topEnd) ,  //Baruch
  //         isRightEnd = curve->is_right_end(topEnd);
  //    if ( !isLeftEnd && !isRightEnd)
  //    {
  //      topEndEvent->add_curve_to_left(curve/*, m_prevPos*/);
  //      topEndEvent->add_curve_to_right(curve);
  //    //  topEndEvent->mark_internal_intersection_point();
  //    }
  //    else
  //      if ( !curve->is_left_end(topEnd)) {
  //        topEndEvent->add_curve_to_left(curve/*, m_prevPos*/);
  //      }
  //      else
  //        if ( ! curve->is_right_end(topEnd)) {
  //          topEndEvent->add_curve_to_right(curve);
  //    }
  //     return true;
  //  } 
  //  
  //  // handle a curve that goes through the bottom point of the vertical curve
  //  const Point_2 &currentPoint = m_currentEvent->get_point();
  //  if (m_traits->point_in_x_range((curve)->get_curve(), currentPoint) &&
  //      m_traits->curve_compare_y_at_x(currentPoint, (curve)->get_curve()) ==
  //      EQUAL)
  //  {
  //     bool isLeftEnd = curve->is_left_end(currentPoint) , 
  //         isRightEnd = curve->is_right_end(currentPoint);
  //     if ( !isLeftEnd && !isRightEnd)
  //     {
  //       m_currentEvent->add_curve_to_left(curve/*, m_prevPos*/);
  //       m_currentEvent->add_curve_to_right(curve);
  //      // m_currentEvent->mark_internal_intersection_point();
  //     }
  //     else
  //       if ( !(curve)->is_left_end(currentPoint)) {
  //        m_currentEvent->add_curve_to_left(curve/*, m_prevPos*/);
  //       }
  //       else
  //         if ( ! (curve)->is_right_end(currentPoint)) {
  //           m_currentEvent->add_curve_to_right(curve);
  //         }
  //    return true;
  //  }
  //  return false;
  //}
  
 /* void handle_vertical_overlap_curves()
  {
    SL_DEBUG(std::cout << "\nhandle_vertical_overlap_curves... ("
                       << m_currentEvent->get_point() << ")";);

    if ( !m_currentEvent->does_contain_vertical_curve() ) {
      SL_DEBUG(std::cout << "no vertical - exiting\n";);
      return;
    }
    SL_DEBUG(std::cout << "\n";);
    SL_DEBUG(PrintVerticals(););

    const Point_2 &point = m_currentEvent->get_point();
    SubCurveListIter iter = m_verticals.begin();
    while ( iter != m_verticals.end() )
    {
      Subcurve *curve = *iter;
      
      if (m_traits->point_in_x_range(curve->get_curve(), point) &&
          m_traits->curve_compare_y_at_x(point, curve->get_curve()) == LARGER)
      {
        iter = m_verticals.erase(iter);
        
      } else if (!curve->is_end_point(point)) {
        EventQueueIter eventIter = m_queue->find(curve->get_top_end());
        CGAL_assertion(eventIter != m_queue->end());
        (eventIter->second)->add_vertical_curve_x_event(m_currentEvent);
        ++iter;
      } else {
        ++iter;
      }
    }

    VerticalCurveList &vcurves = m_currentEvent->get_vertical_curves();
    VerticalCurveListIter vciter = vcurves.begin();
    while ( vciter != vcurves.end() )
    {
      Subcurve *vcurve = *vciter;
      if ( vcurve->is_bottom_end(point) ) {
        m_verticals.push_back(vcurve);
      }
      ++vciter;
    }
  }*/


  // reverse = false ==> we check if the curve starts at the list of curves
  // reverse = true ==> we check if any of the curves in the list start at 
  // the single curve
  void intersect_curve_group(Subcurve *c1, SubCurveList &mylist,
                             PM &pm, bool reverse=false)
  {
    m_tmpOut.clear();
    SL_DEBUG(std::cout << "intersect_curve_group (with out)\n";);
    SL_DEBUG(std::cout << "Intersecting with " << mylist.size()
                       << " curves\n";);
    SubCurveListIter i = mylist.begin();
    SubCurve *prevSubCurve = 0; //the last SubCurve that was handled; 
    bool are_overlap = false; // are current SubCurve and previous SubCurve overlap
   

    while ( i != mylist.end())
    {
     
      X_monotone_curve_2 a , b; 



      //check overlaping of current and previous SubCurve
      if(prevSubCurve && do_curves_overlap(prevSubCurve,*i) )
          are_overlap = true;

      bool flag;
      if ( reverse )
      { 
        flag = CurveStartsAtCurve(*i, c1);
        if ( flag && (c1->get_last_point() != m_currentEvent->get_point()) ) 
        {
          SL_DEBUG(std::cout << "CurveStartsAtCurve 3 \n";);
          m_currentEvent->add_curve_to_right(c1);
          m_currentEvent->add_curve_to_left(c1/*, m_prevPos*/);
         // m_currentEvent->mark_internal_intersection_point(); //Baruch
          SL_DEBUG(std::cout << "splitting " << (c1)->get_last_curve()
                             << " at " 
                             << m_currentEvent->get_point() << "\n";);
          if ( (c1)->is_source_left_to_target() ) 
            m_traits->curve_split((c1)->get_last_curve(), a, b, 
                                  m_currentEvent->get_point());
          else
            m_traits->curve_split((c1)->get_last_curve(), b, a, 
                                  m_currentEvent->get_point());

          if(!are_overlap)
          {
            Halfedge_handle h(NULL);
            h = insert_to_pm(a, c1, h, pm);
          }
          else
            if(c1->get_last_event()->get_insert_info()->dec_right_curves_counter() == 0)
            {
              //delete c1->get_last_event();
                m_eventAlloc.destroy(c1->get_last_event());
                m_eventAlloc.deallocate(c1->get_last_event(),1);
            }

          (c1)->set_last_point(m_currentEvent->get_point());
          (c1)->set_last_curve(b); 
          //(c1)->set_last_subcurve(a); 
          (c1)->set_last_event(m_currentEvent);
          m_tmpOut.push_back(c1);
        }
      }
      else   // !reverse
      {
        /*X_monotone_curve_2& a = (*i)->get_last_subcurve_ref();
        X_monotone_curve_2& b = (*i)->get_last_curve_ref();*/
        flag = CurveStartsAtCurve(c1, *i);
        if ( flag && ((*i)->get_last_point() != m_currentEvent->get_point())) 
        {
          SL_DEBUG(std::cout << "CurveStartsAtCurve 3 \n";);
          m_currentEvent->add_curve_to_right(*i);
          m_currentEvent->add_curve_to_left(*i/*, m_prevPos*/);
        //  m_currentEvent->mark_internal_intersection_point(); //Baruch
          SL_DEBUG(std::cout << "splitting " << (*i)->get_last_curve() 
                             << " at " 
                             << m_currentEvent->get_point() << "\n";);
          if ( (*i)->is_source_left_to_target() ) 
            m_traits->curve_split((*i)->get_last_curve(), a, b, 
                                  m_currentEvent->get_point());
          else
            m_traits->curve_split((*i)->get_last_curve(), b, a, 
                                  m_currentEvent->get_point());

          if(!are_overlap)
          {
            Halfedge_handle h(NULL);
            h = insert_to_pm(a, *i, h, pm);
          }
          else
            if((*i)->get_last_event()->get_insert_info()->dec_right_curves_counter() == 0)
            {
              //delete (*i)->get_last_event();
               m_eventAlloc.destroy((*i)->get_last_event());
               m_eventAlloc.deallocate((*i)->get_last_event(),1);
            }

          (*i)->set_last_point(m_currentEvent->get_point());
          (*i)->set_last_curve(b); 
          //(*i)->set_last_subcurve(a); 
          (*i)->set_last_event(m_currentEvent);
          m_tmpOut.push_back(*i);
          
        }
      }
      
      intersect(c1, *i);
      prevSubCurve= *i; //update current SubCurve to be the previous
      ++i;
       are_overlap = false; 
    }    
  }

private:

  /*! Loop over the curves to the right of the status line and handle them:
   * - if we are at the beginning of the curve, we insert it to the status 
   *   line, then we look if it intersects any of its neighbours.
   * - if we are at an intersection point between two curves, we add them
   *   to the status line and attempt to intersect them with their neighbours. 
   * - We also check to see if the two intersect again to the right of the point.
   */
  void handle_right_curves(PM &pm, SweepLinePlanarmap &tag)
  {
    SL_DEBUG(std::cout << "Handling right curves (" ;);
    SL_DEBUG(std::cout << m_currentEvent->get_point() << ")\n";);
    int numRightCurves = m_currentEvent->get_num_right_curves();
    if ( numRightCurves == 0 ) //no right cures, return
      return;
    
    //m_currentPos = m_sweepLinePos;
    if ( numRightCurves == 1 )
    {
      SL_DEBUG(std::cout << " - beginning of curve " << std::endl;);
        
      SL_DEBUG(
               Subcurve *tmp1 = *(m_currentEvent->right_curves_begin());
               PRINT_INSERT(tmp1);
               );
        
      StatusLineIter slIter = 
        m_statusLine->insert(m_status_line_insert_hint, 
                             *(m_currentEvent->right_curves_begin()));
      
      (*(m_currentEvent->right_curves_begin()))->set_hint(slIter); //xxx
      m_status_line_insert_hint = slIter; ++m_status_line_insert_hint;
      
      SL_DEBUG(PrintStatusLine(););
        
      // if this is the only curve on the status line, nothing else to do
      if ( m_statusLine->size() == 1 )
        return; 
  
      StatusLineIter prev = slIter; // the previous neighbour of the curve at the status line
      StatusLineIter next = slIter; // the next neighbour of the curve at the status line
      ++next;
      
      // 'mylist' will hold the two neighbours of the curve ,
      //  and all of their overlapped curves
      SubCurveList mylist;
      if ( slIter != m_statusLine->begin() )
      {
        --prev;
        StatusLineIter tmp = prev;
        mylist.push_back(*prev); // push the previous neighbour 

        //find all of the curves that overlap with *prev and push them to 'mylist'
        while ( tmp != m_statusLine->begin() ) 
        {
          --tmp;
          if ( do_curves_overlap(*prev, *tmp) )
            mylist.push_back(*tmp);
          else 
            break;
        }
      }
      
      if ( next != m_statusLine->end() )
      { 
        StatusLineIter tmp = next;
        mylist.push_back(*next); // push the next neighbour
        ++tmp;

        //find all of the curves that overlap with *next and push them to 'mylist'
        while ( tmp != m_statusLine->end() ) 
        {
          if ( do_curves_overlap(*next, *tmp) )
          {
            mylist.push_back(*tmp);
            ++tmp;
          }
          else
            break;
        }
      }
      intersect_curve_group(*(m_currentEvent->right_curves_begin()), mylist, pm);
      
    } 
    else  // if we've reached here , numRightCurves > 1 
    {
      /* this block takes care of 
      //
      //           /
      //          /
      //       --------
      //          \
      //           \
      */
      int numLeftCurves = m_currentEvent->get_num_left_curves();
      if ( numLeftCurves == 0 ) 
      {

        SL_DEBUG(std::cout << " - handling special case " << std::endl;);

        StatusLineIter slIter;
        EventCurveIter currentOne = m_currentEvent->right_curves_begin();
        while ( currentOne != m_currentEvent->right_curves_end() ) 
        {
          //find an iterator to the first Subcurve in the statusLine
          // that is equal to or greater than currentOne

          slIter = m_statusLine->lower_bound(*currentOne);
          if ( slIter != m_statusLine->end() )
          {
            Subcurve *c = *slIter;
            if ( CurveStartsAtCurve(*currentOne, c)) 
            {
              m_currentEvent->add_curve_to_left(c/*, m_sweepLinePos*/); 
              m_currentEvent->add_curve_to_right(c);
            //  m_currentEvent->mark_internal_intersection_point(); //Baruch
              X_monotone_curve_2 a, b;
              if ( c->is_source_left_to_target() ) 
              {
                m_traits->curve_split(c->get_last_curve(), a, b, 
                                      m_currentEvent->get_point());
              } 
              else
              {
                m_traits->curve_split(c->get_last_curve(), b, a, 
                                      m_currentEvent->get_point());
              }
              
              Halfedge_handle h(NULL);
              insert_to_pm(a, c, h, pm);
              
              c->set_last_point(m_currentEvent->get_point());
              c->set_last_curve(b);
              c->set_last_event(m_currentEvent);
              break;
            }
          }
          currentOne++;
        }
      } // end block ... ( numLeftCurves == 0 )
      
      SubCurveList mylist;
      SubCurveList prevlist;
      SubCurveList currentlist;
      
      SL_DEBUG(std::cout << " - intersection point " << std::endl;);
      EventCurveIter firstOne = m_currentEvent->right_curves_begin();
      EventCurveIter lastOne = m_currentEvent->right_curves_end(); --lastOne;
      EventCurveIter rightCurveEnd = m_currentEvent->right_curves_end();
      
      PRINT_INSERT(*firstOne);
      
      StatusLineIter slIter = m_statusLine->insert(m_status_line_insert_hint, 
                                                   *firstOne);
      (*firstOne)->set_hint(slIter);//xxx
      
      SL_DEBUG(PrintStatusLine(););
      if ( slIter != m_statusLine->begin() )
      { 
        StatusLineIter prev = slIter; --prev;
          
        // find all curves that are overlapping with the prev curve
        StatusLineIter tmp = prev;
        prevlist.push_back(*prev);
        while ( tmp != m_statusLine->begin() ) 
        {
          --tmp;
          if ( do_curves_overlap(*prev, *tmp))
            prevlist.push_back(*tmp);
          else
            break;
        }
 
        intersect_curve_group(*slIter, prevlist, pm);
      }
      currentlist.push_back(*firstOne);

      EventCurveIter currentOne = firstOne; ++currentOne;
      EventCurveIter prevOne = firstOne;

      while ( currentOne != rightCurveEnd )
      {
        //m_currentPos = m_sweepLinePos;
        PRINT_INSERT(*currentOne);
       // ++slIter;
      //  slIter = m_statusLine->insert(slIter, *currentOne);
        slIter = m_statusLine->insert(m_status_line_insert_hint, *currentOne); // YYY
        (*currentOne)->set_hint(slIter);//xxx
        
        SL_DEBUG(PrintStatusLine(););
        if ( do_curves_overlap(*currentOne, *prevOne))
        {
          intersect_curve_group(*currentOne, currentlist, pm);
          currentlist.push_back(*currentOne);
        }
        else
        {
          prevlist = currentlist;
          currentlist.clear();
          currentlist.push_back(*currentOne);
        }
        
        intersect_curve_group(*currentOne, prevlist, pm);
        prevOne = currentOne;
        ++currentOne;
      }
     // m_status_line_insert_hint = slIter; ++m_status_line_insert_hint;
    
      lastOne = currentOne; --lastOne;
     // m_currentPos = m_sweepLinePos;
      PRINT_INSERT(*lastOne);
    
      SL_DEBUG(PrintStatusLine(););
      StatusLineIter next = slIter; ++next;
      if ( next != m_statusLine->end() ) {
        intersect_curve_group(*next, currentlist, pm, true);
        StatusLineIter tmp = next; ++tmp;
        while ( tmp != m_statusLine->end() ) 
        {
          if ( do_curves_overlap(*next, *tmp))
          {
            intersect_curve_group(*tmp, currentlist, pm, true);
            ++tmp;
          }
          else
            break;
        }
      }
    }
  }


  /*! Insert a curve to the planar map.
   *  If an identical curve was already inserted into the planarmap, it is 
   *  not inserted again. 
   *
   *  @param cv the curve to insert
   *  @param leftCurve the original curve
   *  @param hhandle a prev halfedge handle (may be NULL)
   *  @param pm a reference to the planar map
   */
  Halfedge_handle insert_to_pm(const X_monotone_curve_2 &cv, 
                               SubCurve *leftCurve,
                               Halfedge_handle hhandle, PM &pm)
  {
    SL_DEBUG(std::cout << "*X inserting " << cv << "(" 
                       << leftCurve << ")\n";);
    static SubCurve *prevCurve = 0;
    static X_monotone_curve_2 prevXCv;
    
    Event *lastEvent = leftCurve->get_last_event();
    PmwxInsertInfo *insertInfo = lastEvent->get_insert_info();
    
    SL_DEBUG(std::cout << "lastEvent = " << lastEvent << "\n";
             lastEvent->Print();
             insertInfo->Print(););
    
    // if this is the same as the previous curve, don't add it again
    if ( prevCurve && similar_curves(cv, prevXCv)) {
      leftCurve->set_last_event(m_currentEvent);

      if ( m_change_not ) {
        m_change_not->add_edge(cv, hhandle, true, true);
      }
      if(lastEvent->get_insert_info()->dec_right_curves_counter() == 0)
      {
        //delete lastEvent;
         m_eventAlloc.destroy(lastEvent);
         m_eventAlloc.deallocate(lastEvent,1);
      }

      return hhandle;
    }
    prevCurve = leftCurve;
    prevXCv = cv;
    
    Halfedge_handle res; 
    SL_DEBUG(std::cout << "get_halfedge_jump_count : curve is " 
                       << leftCurve->get_curve() << "\n";);
    int jump = lastEvent->get_halfedge_jump_count(leftCurve);

    Point_2 p1, p2;
    if (m_traits->compare_xy(m_traits->curve_source(cv), m_traits->curve_target(cv)) == SMALLER) {
      p1 = m_traits->curve_source(cv);
      p2 = m_traits->curve_target(cv);
    } else {
      p2 = m_traits->curve_source(cv);
      p1 = m_traits->curve_target(cv);
    }
      
    // if the previous event on the curve is not in the planar map yet
    if ( insertInfo->get_halfedge_handle() == Halfedge_handle(NULL) ) 
    {
      // we have a handle from the previous insert
      if ( hhandle != Halfedge_handle(NULL) )
      {
        if ( !m_traits->point_equal( hhandle->target()->point(), p2 ))
          hhandle = hhandle->twin();
        SL_DEBUG(std::cout << "  from vertex (1)";
                 std::cout << hhandle->source()->point() << " " 
                 << hhandle->target()->point() << "\n";);
        
                
       res = pm.non_intersecting_insert_from_vertex(cv, hhandle, m_change_not);
                        

      }
     
      else
      {
        // if this is the first left curve being inserted
        SL_DEBUG(std::cout << "  in face interior\n";);
        res = pm.insert_in_face_interior(cv, pm.unbounded_face(), m_change_not);
      }
    } else 
      // the previous event on the curve is already in the planar map. 
      // Let's use it.
    {
      Halfedge_handle prev = insertInfo->get_halfedge_handle();
      if ( !m_traits->point_equal( prev->target()->point(), p1 ))
        prev = prev->twin();
      // skip to the right halfedge
      SL_DEBUG(std::cout << "Skipping " << jump << " steps\n";);
      for ( int i = 0 ; i < jump ; i++ )
        prev = (prev->next_halfedge())->twin();
      
      // we have a handle from the previous insert
      if ( hhandle != Halfedge_handle(NULL) ) 
      {
        if ( !m_traits->point_equal( hhandle->target()->point(), p2 ))
          hhandle = hhandle->twin();

        CGAL_assertion(prev->face() == hhandle->face());
        SL_DEBUG(std::cout << "  at vertices";
                 std::cout << prev->source()->point() << " " 
                           << prev->target()->point() << " --- ";
                 std::cout << hhandle->source()->point() << " " 
                           << hhandle->target()->point() << "\n";);
         
          res = pm.non_intersecting_insert_at_vertices(cv, prev, hhandle, 
               m_change_not);

      } else {

        SL_DEBUG(std::cout << "  from vertex (2)";
                 std::cout << prev->source()->point() << " " 
                           << prev->target()->point() << "\n";);

        res = pm.non_intersecting_insert_from_vertex(cv, prev, m_change_not); 
                                
      }
    }
  
    SL_DEBUG(std::cout << "*** returning: (" << res->source()->point() << " " 
                       << res->target()->point()  << ")\n\n";);
   
   
    // update the information in the events so they can be used in the future
    /*bool exist_vertical = insertInfo->get_vertical_below_event_flag() ||
                          insertInfo->get_vertical_above_event_flag();*/

    if ( lastEvent->get_num_left_curves() == 0 &&  
         lastEvent->is_curve_largest(leftCurve) /*&& !exist_vertical */)
    {
       insertInfo->set_halfedge_handle(res->twin());
    }
    insertInfo = m_currentEvent->get_insert_info();
    insertInfo->set_halfedge_handle(res);

    if(lastEvent->get_insert_info()->dec_right_curves_counter() == 0)
    {
      //m_events.remove(lastEvent);
      //delete lastEvent;
       m_eventAlloc.destroy(lastEvent);
       m_eventAlloc.deallocate(lastEvent,1);
    }
    return res;

  }

 /* void insert_to_pm_v(const X_monotone_curve_2 &a, SubCurve *origCurve, 
                      Event *topEvent, Event *bottomEvent, PM &pm);*/


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
//template <class CurveInputIterator, class SweepLineTraits_2, 
//          class PM_, class Change_notification_>
//void 
//Pmwx_aggregate_insert<CurveInputIterator, SweepLineTraits_2,
//                      PM_, Change_notification_>::
//insert_to_pm_v(const X_monotone_curve_2 &cv, SubCurve *origCurve, 
//               Event *topEvent, Event *bottomEvent, PM &pm)
//{
//  SL_DEBUG(std::cout << "insert_to_pm_v \n";
//           std::cout << "*V inserting " << cv << ")\n";);
//
//  PmwxInsertInfo *topII = topEvent->get_insert_info();
//  PmwxInsertInfo *bottomII = bottomEvent->get_insert_info();
//  Halfedge_handle res; 
//
//  // if the curve is already in the planar map, update the data and return
//  if (vertical_subcurve_exists(cv)) {
//    origCurve->set_last_event(m_currentEvent);
//    SL_DEBUG(std::cout << "*** returning (curve already exists\n\n";);
//    return;
//  }
//
//  m_verticalSubCurves.push_back(cv);
//
//  if ( topII->get_halfedge_handle() == Halfedge_handle(NULL))
//  {
//    if ( bottomII->get_halfedge_handle() == Halfedge_handle(NULL))
//    {
//      SL_DEBUG(std::cout << "  in face interior\n";);
//      
//      res = pm.insert_in_face_interior(cv, pm.unbounded_face(), m_change_not);
//      if ( !origCurve->is_source_left_to_target() ){
//        res = res->twin();
//      }        
//    } else 
//    {
//      SL_DEBUG(std::cout << "  from vertex (1) ";
//               std::cout << bottomII->get_halfedge_handle()->source()->point() 
//                         << bottomII->get_halfedge_handle()->target()->point()
//                         << "\n";);
//      Halfedge_handle h1 = bottomII->get_halfedge_handle();
//      if ( h1->target()->point() != m_traits->curve_source(cv) &&
//           h1->target()->point() != m_traits->curve_target(cv) )
//        h1 = h1->twin();
//     res = pm.non_intersecting_insert_from_vertex(cv, h1, m_change_not);  
//                        
//    }
//  } else 
//  {
//    if ( bottomII->get_halfedge_handle() == Halfedge_handle(NULL))
//    {
//      SL_DEBUG(std::cout << "  from vertex (2) ";
//               std::cout << topII->get_halfedge_handle()->source()->point() 
//                         << topII->get_halfedge_handle()->target()->point() 
//                         << "\n";);
//      Halfedge_handle h1 = topII->get_halfedge_handle();
//      if ( h1->target()->point() != m_traits->curve_source(cv) &&
//           h1->target()->point() != m_traits->curve_target(cv) )
//        h1 = h1->twin();
//      res = pm.non_intersecting_insert_from_vertex(cv, h1, m_change_not);
//                        
//      res = res->twin();
//    } else 
//    {
//      SL_DEBUG(std::cout << "  at vertices";
//               std::cout << bottomII->get_halfedge_handle()->source()->point();
//               std::cout << bottomII->get_halfedge_handle()->target()->point();
//               std::cout << topII->get_halfedge_handle()->source()->point();
//               std::cout << topII->get_halfedge_handle()->target()->point() 
//                         << "\n";);
//      Halfedge_handle h1 = bottomII->get_halfedge_handle();
//      Halfedge_handle h2 = topII->get_halfedge_handle();
//      if ( h1->target()->point() != m_traits->curve_source(cv) &&
//           h1->target()->point() != m_traits->curve_target(cv) )
//        h1 = h1->twin();
//      if ( h2->target()->point() != m_traits->curve_source(cv) &&
//           h2->target()->point() != m_traits->curve_target(cv) )
//        h2 = h2->twin();
//      
//      CGAL_assertion(h1->face() == h2->face());
//      res = pm.non_intersecting_insert_at_vertices(cv, h1, h2, m_change_not);
//                        
//    }
//  }
//
//   
//  if ( topEvent->get_num_left_curves() == 0  ) 
//  {
//    if(topII->get_vertical_above_event_flag())/*res->twin() != res->next_halfedge()*/  
//    {
//       topII->set_halfedge_handle(res->next_halfedge()->twin()); 
//       topII->set_vertical_above_event_flag(true);
//    }
//    else  
//    {
//      topII->set_halfedge_handle(res);  
//    }
//  } 
//  topII->set_vertical_below_event_flag(true);
//  bottomII->set_halfedge_handle(res->twin());
//  bottomII->set_vertical_above_event_flag(true);
//
//  /* if(bottomII->dec_right_curves_counter() == 0)
//    {
//      std::cout<< "deleting event bottom ...\n";
//      std::cout.flush();
//      delete bottomEvent;
//      std::cout<< "finished deleting event bottom ...\n";
//      std::cout.flush();
//    }*/
//
//  SL_DEBUG(std::cout << "*** returning: (" << res->source()->point() << " " 
//                     << res->target()->point()  << ")\n\n";);
//
//  
//}
  
CGAL_END_NAMESPACE

#endif // CGAL_PMWX_AGGREGATE_INSERT_H
