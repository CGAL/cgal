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
  typedef typename PM::Halfedge_handle    Halfedge_handle;
  typedef typename PM::Edge_iterator      Edge_iterator;

  

  typedef Pmwx_sweep_line_curve<SweepLineTraits_2, 
                                Halfedge_handle> SubCurve;
  typedef Pmwx_sweep_line_event<SweepLineTraits_2, SubCurve> Event;
  typedef typename SubCurve::PmwxInsertInfo PmwxInsertInfo;

  typedef Change_notification_ Change_notification;
  typedef Sweep_line_2_impl<CurveInputIterator, Traits, Event, SubCurve> Base;

  typedef Pmwx_sweep_line_curve<SweepLineTraits_2, 
                                typename PM_::Halfedge_handle>  Subcurve;

  // repeated typedefs from the base class to avoid warnings
  typedef typename Base::EventQueueIter EventQueueIter;
  typedef typename Base::EventCurveIter EventCurveIter;
  typedef typename Base::StatusLineIter StatusLineIter;
  typedef typename Base::SubCurveList SubCurveList;
  typedef typename Base::SubCurveListIter SubCurveListIter;
  typedef typename Base::SweepLinePlanarmap SweepLinePlanarmap;
  typedef typename Base::EventQueueValueType EventQueueValueType;

#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_3
  using Base::m_queue;
  using Base::m_traits;
  using Base::m_sweepLinePos;
  using Base::m_currentEvent;
  using Base::m_statusLine;
  using Base::m_tmpOut;
  using Base::m_status_line_insert_hint;
#endif

  Pmwx_aggregate_insert() : 
    Base(), m_change_not(NULL)
    {}
  
  Pmwx_aggregate_insert(Traits *traits_) : 
    Base(traits_), m_change_not(NULL) 
    {} 
  
  virtual ~Pmwx_aggregate_insert() {}
  

  /*! Initializes the data structures to work with:
    - x-monotonize the input curves
    - for each end point of each curve create an event
    - for each curve in the planarmap do the same
    - initialize the event queue
    -
  */  
  void init(CurveInputIterator begin, CurveInputIterator end, PM &pm)
  {
    Edge_iterator eit;

    for (eit = pm.edges_begin(); eit != pm.edges_end(); ++eit) 
    {
      m_xcurves.push_back(eit->curve());
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

    while ( eventIter != m_queue->end() )
    {
      const Point_2 &p = (*eventIter).first;
      m_sweepLinePos = p;
      m_currentEvent = eventIter->second;
      SL_DEBUG(std::cout << "------------- " << p << " --------------"
                         << std::endl;
               PrintStatusLine();
               m_currentEvent->Print();
               );

      handle_left_curves(pm, tag);

      m_queue->erase(eventIter);
      handle_right_curves(pm, tag);
      if(m_currentEvent->get_insert_info()->get_right_curves_counter() == 0)
      {
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
    const Point_2 &eventPoint = m_currentEvent->get_point();

    Halfedge_handle h(NULL);
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
      
        leftCurvePrev = *leftCurveIter; 
        ++leftCurveIter;
        are_overlap = false;
        
    }  
    
    // when done handling the left curves, we prepare for the right curves
   // m_currentEvent->init_right_curves();
  }
  /////////////////////////////////////////////////////////////////////////////
  void handle_left_curves_no_overlap(PM &pm, SweepLinePlanarmap &tag)
  {
    EventCurveIter leftCurveIter = m_currentEvent->left_curves_begin();
    const Point_2 &eventPoint = m_currentEvent->get_point();

    Halfedge_handle h(NULL);
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
      
        ++leftCurveIter;  
    }  
    
    // when done handling the left curves, we prepare for the right curves
   // m_currentEvent->init_right_curves();
  }
  ////////////////////////////////////////////////////////////////////////


 
 
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
          m_currentEvent->add_curve_to_left(c1);
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
                m_eventAlloc.destroy(c1->get_last_event());
                m_eventAlloc.deallocate(c1->get_last_event(),1);
            }

          (c1)->set_last_point(m_currentEvent->get_point());
          (c1)->set_last_curve(b); 
          (c1)->set_last_event(m_currentEvent);
          m_tmpOut.push_back(c1);
        }
      }
      else   // !reverse
      {
        flag = CurveStartsAtCurve(c1, *i);
        if ( flag && ((*i)->get_last_point() != m_currentEvent->get_point())) 
        {
          SL_DEBUG(std::cout << "CurveStartsAtCurve 3 \n";);
          m_currentEvent->add_curve_to_right(*i);
          m_currentEvent->add_curve_to_left(*i);
   
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
               m_eventAlloc.destroy((*i)->get_last_event());
               m_eventAlloc.deallocate((*i)->get_last_event(),1);
            }

          (*i)->set_last_point(m_currentEvent->get_point());
          (*i)->set_last_curve(b);  
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
              m_currentEvent->add_curve_to_left(c); 
              m_currentEvent->add_curve_to_right(c);
      
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
        PRINT_INSERT(*currentOne);
        slIter = m_statusLine->insert(m_status_line_insert_hint, *currentOne); 
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
    
      lastOne = currentOne; --lastOne;
    
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
   
   
    if ( lastEvent->get_num_left_curves() == 0 &&  
         lastEvent->is_curve_largest(leftCurve) )
    {
       insertInfo->set_halfedge_handle(res->twin());
    }
    insertInfo = m_currentEvent->get_insert_info();
    insertInfo->set_halfedge_handle(res);

    if(lastEvent->get_insert_info()->dec_right_curves_counter() == 0)
    {
       m_eventAlloc.destroy(lastEvent);
       m_eventAlloc.deallocate(lastEvent,1);
    }
    return res;

  }



  Change_notification *m_change_not;
};


  
CGAL_END_NAMESPACE

#endif // CGAL_PMWX_AGGREGATE_INSERT_H
