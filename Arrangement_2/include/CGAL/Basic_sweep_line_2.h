// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// $URL$
// $Id$
// 
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 (based on old version by Tali Zvi)

#ifndef CGAL_BASIC_SWEEP_LINE_2_H
#define CGAL_BASIC_SWEEP_LINE_2_H

#include <vector>
#include <algorithm>
#include <iterator>
#include <CGAL/assertions.h>
#include <CGAL/memory.h>
#include <CGAL/Sweep_line_2/Sweep_line_functors.h>
#include <CGAL/Sweep_line_2/Sweep_line_subcurve.h>
#include <CGAL/Sweep_line_2/Sweep_line_event.h>
#include <CGAL/Multiset.h>
#include <CGAL/Arrangement_2/Arr_traits_adaptor_2.h>

#ifndef CGAL_DEBUG_SWEEP_LINE

#define CGAL_SL_DEBUG(a)
#define CGAL_PRINT_INSERT(a)
#define CGAL_PRINT_ERASE(a)
#define CGAL_PRINT_NEW_EVENT(p, e) 
#define CGAL_PRINT(a)


#else

#define CGAL_SL_DEBUG(a) {a}
#define CGAL_PRINT_INSERT(a) { std::cout << "+++ inserting "; \
                          (a)->Print(); \
                          std::cout << "    currentPos = "; \
                          PrintEvent(this->m_currentEvent); \
                          std::cout << "\n"; \
                          }
#define CGAL_PRINT_ERASE(a)  { std::cout << "--- erasing " ; \
                          (a)->Print(); }
#define CGAL_PRINT_NEW_EVENT(p, e) \
{ std::cout << "%%% a new event was created at " << (p) << std::endl; \
  (e)->Print(); }
#define CGAL_PRINT(a) { std::cout << a ; }

#endif



CGAL_BEGIN_NAMESPACE

/*!
  Basic_Sweep_line_2 is a class that implements the sweep line algorithm
  for non-intersecting x-monotone curves.
  It extends the algorithm to support not only segments but general curves
  as well and isolated points.
  The curves are defined by the traits class that is one of the template 
  arguments.

  The algorithm is also extended to support the following degenerate cases:
  - vertical segments

  General flow:
  After the initialization stage, the events are handled from left to right.

  For each event

    Handle left curves - iterate over the curves that intersect 
                 at the event point and defined to the left of the 
                 event. 
    Handle right curves - iterate over the curves that intersect 
                 the event point and defined to the right of the 
                 event point.
  End

  Convensions through out the code:
  In order to make the code as readable as possible, some convensions were 
  made in regards to variable naming:

    slIter - an iterator to the status line, always points to a curve.

*/
template < class Traits_,
           class SweepVisitor,
           class CurveWrap = Sweep_line_subcurve<Traits_>,
           typename SweepEvent = Sweep_line_event<Traits_, CurveWrap>,
           typename Allocator = CGAL_ALLOCATOR(int) >
class Basic_sweep_line_2
{
public:

  typedef Traits_                                        Traits;
  typedef Arr_traits_basic_adaptor_2<Traits_>            Traits_adaptor;
  typedef typename Traits_adaptor::Point_2                       Point_2;
  typedef typename Traits_adaptor::X_monotone_curve_2            X_monotone_curve_2;

  typedef SweepEvent                                     Event;
  typedef Event_less_functor<Traits_adaptor, Event>              EventLess;
  typedef Multiset<Event*, EventLess, Allocator>         EventQueue; 
  typedef typename EventQueue::iterator                  EventQueueIter;

  typedef typename Event::SubCurveIter                   EventCurveIter;

  typedef Sweep_line_event<Traits, CurveWrap>            Base_event;
  typedef typename Base_event::Attribute                 Attribute;
  
  typedef CurveWrap                                      Subcurve;

  typedef Sweep_line_subcurve<Traits>                    Base_subcurve;
  typedef Status_line_curve_less_functor<Traits,
                                         Base_subcurve>  StatusLineCurveLess;
  typedef Multiset<Base_subcurve*,
                   StatusLineCurveLess,
                   Allocator>                            StatusLine;
  typedef typename StatusLine::iterator                  StatusLineIter;

  typedef typename Allocator::template rebind<Event>     EventAlloc_rebind;
  typedef typename EventAlloc_rebind::other              EventAlloc;

  typedef typename Allocator::template rebind<Subcurve>  SubcurveAlloc_rebind;
  typedef typename SubcurveAlloc_rebind::other           SubCurveAlloc;

protected:

  struct CompEventPtr
  {
    Comparison_result operator() (Event *e1, Event *e2) const
    {
      if (e1 < e2)
        return (SMALLER);
      if (e1 > e2)
        return (LARGER);
      return (EQUAL);
    }
  };

  typedef Multiset<Event*, CompEventPtr>           Allocated_events_set;
  typedef typename Allocated_events_set::iterator  Allocated_events_iterator;

public:

  /*!
   * Constructor.
   * \param visitor A pointer to a sweep-line visitor object.
   */
  Basic_sweep_line_2 (SweepVisitor* visitor) :
      m_traits(new Traits_adaptor()),
      m_traitsOwner(true),
      m_statusLineCurveLess(m_traits, &m_currentEvent),
      m_queueEventLess(m_traits),
      m_queue(new EventQueue(m_queueEventLess)),
      m_statusLine(m_statusLineCurveLess),
      m_status_line_insert_hint(m_statusLine.begin()),
      m_num_of_subCurves(0),
      m_visitor(visitor)
  {
    m_visitor->attach(this);
  }


  /*!
   * Constructor.
   * \param traits A pointer to a sweep-line traits object.
   * \param visitor A pointer to a sweep-line visitor object.
   */
  Basic_sweep_line_2(Traits *traits, SweepVisitor* visitor) :
      m_traits(static_cast<Traits_adaptor*>(traits)),
      m_traitsOwner(false),
      m_statusLineCurveLess(m_traits, &m_currentEvent),
      m_queueEventLess(m_traits),
      m_queue(new EventQueue(m_queueEventLess)),
      m_statusLine(m_statusLineCurveLess),
      m_status_line_insert_hint(m_statusLine.begin()),
      m_num_of_subCurves(0),
      m_visitor(visitor)
  {
    m_visitor->attach(this);
  }

  /*! Destrcutor. */
  virtual ~Basic_sweep_line_2()
  {
    if(m_traitsOwner)
      delete m_traits;
    delete m_queue;

    // Free all the event that have not been de-allocated so far.
    Allocated_events_iterator      iter;
    Event                         *p_event;

    for (iter = m_allocated_events.begin();
         iter != m_allocated_events.end(); ++iter)
    {
      p_event = *iter;
      m_eventAlloc.destroy(p_event);
      m_eventAlloc.deallocate(p_event,1);
    }
  }

  /*!
   * Run the sweep-line with a range of curves.
   * \param curves_begin An iterator for the first curve in the range.
   * \param curves_end A past-the-end iterator for the range.
   * \pre The value-type of CurveInputIterator is X_monotone_curve_2.
   */
  template<class CurveInputIterator>
  void sweep(CurveInputIterator curves_begin,
             CurveInputIterator curves_end)
  {
    m_visitor->before_sweep();
    _init_sweep(curves_begin, curves_end);
    m_visitor ->after_init();
    _sweep();
    _complete_sweep();
     m_visitor ->after_sweep();
  }


   /*!
    * Run the sweep-line with a range of x-monotone curves and a range   
    * of action event points(if a curve passed through an action point,it will
    * be splitted).
    * \param curves_begin  An iterator for the first x-monotone curve in the
    *                      range.
    * \param curves_end A past-the-end iterator for this range.
    * \param points_begin An iterator for the first point in the range.
    * \param points_end A past-the-end iterator for this range.
    * \pre The value-type of XCurveInputIterator is the traits' 
    *      X_monotone_curve_2, and the value-type of PointInputIterator is the
    *      traits' Point_2.
    */
  template<class CurveInputIterator, class PointInputIterator>
  void sweep (CurveInputIterator curves_begin,
              CurveInputIterator curves_end,
	            PointInputIterator action_points_begin,
              PointInputIterator action_points_end)
  {
    m_visitor->before_sweep();
    _init_sweep(curves_begin, curves_end);
    _init_points(action_points_begin, action_points_end, Base_event::ACTION);
    m_visitor ->after_init();
    _sweep();
    _complete_sweep();
     m_visitor ->after_sweep();
  }

   /*!
    * Run the sweep-line with a range of x-monotone curves, a range   
    * of action event points (if a curve passed through an action point,it will
    * be splitted) and a range of query points (if a curve passed through a
    * query point,it will not be splitted).
    * \param curves_begin  An iterator for the first x-monotone curve in the
    *                      range.
    * \param curves_end A past-the-end iterator for this range.
    * \param points_begin An iterator for the first point in the range.
    * \param points_end A past-the-end iterator for this range.
    * \pre The value-type of XCurveInputIterator is the traits' 
    *      X_monotone_curve_2, and the value-type of PointInputIterator is the 
    *      traits' Point_2.
    */
  template<class CurveInputIterator, class ActionPointItr,class QueryPointItr>
  void sweep (CurveInputIterator curves_begin,
              CurveInputIterator curves_end,
	            ActionPointItr action_points_begin,
              ActionPointItr action_points_end,
              QueryPointItr query_points_begin,
              QueryPointItr query_points_end)
  {
    m_visitor->before_sweep();
    _init_sweep(curves_begin, curves_end);
    _init_points(action_points_begin, action_points_end, Base_event::ACTION);
    _init_points(query_points_begin, query_points_end, Base_event::QUERY);
    m_visitor ->after_init();
    _sweep();
    _complete_sweep();
    m_visitor ->after_sweep();
  }


  /*! Get an iterator for the first subcurve in the status line. */
  StatusLineIter  status_line_begin()
  {
    return m_statusLine.begin();
  }

  /*! Get a past-the-end iterator for the subcurves in the status line. */
  StatusLineIter  status_line_end()
  {
    return m_statusLine.end();
  }

  /*! Get the status line size */
  unsigned int status_line_size() const
  {
    return m_statusLine.size();
  }

  /*! return bool iff m_statusLine is empty */
  bool is_status_line_empty() const
  {
    return (!m_statusLine.size());
  }

  /*! Stop the sweep by erasing the X-strucure (except for the current event)
   *  can be called by the visitor during 'arter_handle_event'.
   */
  void stop_sweep()
  {
    EventQueueIter qiter= this ->m_queue->begin();
    ++qiter;
    for(;
        qiter != this ->m_queue->end();
        ++qiter)
    {
      this ->deallocate_event(*qiter);
    }
    this -> m_statusLine.clear();
    m_status_line_insert_hint = this -> m_statusLine.begin();
  
    CGAL_assertion(!m_queue->empty());
    EventQueueIter second = m_queue->begin(); ++second;
    while(second != m_queue->end())
    {
      EventQueueIter next = second;
      ++next;
      m_queue->erase(second);
      second = next;
    }
    
  }


   /*! Deallocate event object, it is a public method to allow the visitor
    *  to manage the events deallocation (if he wants to) 
    */
  void deallocate_event(Event* event)
  {
    m_allocated_events.erase (event);

    m_eventAlloc.destroy(event);
    m_eventAlloc.deallocate(event,1);
  }

  /*! Get the current event */
  Event* current_event()
  {
    return m_currentEvent;
  }

  /*! Get the traits object */
  Traits* traits()
  {
    return m_traits;
  }

 

  protected:

    /*! Preform the main sweep-line loop. */
  void _sweep()
  {
    // Looping over the events in the queue.
    EventQueueIter eventIter = m_queue->begin();

    while (eventIter != m_queue->end())
    {
      // Get the next event from the queue.
      m_currentEvent = *eventIter;

      CGAL_PRINT("------------- ");
      CGAL_SL_DEBUG(PrintEvent(m_currentEvent););
      CGAL_PRINT ( " --------------\n");
      CGAL_SL_DEBUG(PrintStatusLine();
               m_currentEvent->Print(););
      
      // Handle the subcurves that are to the left of the event point (i.e., 
      // subcurves that we are done with).
      _handle_left_curves();

      // Handle the subcurves to the right of the event point, reorder them
      // and test for intersections between them and their immediate neighbors
      // on the status line.
      _handle_right_curves();

      // Inform the visitor about the event.
      if (m_visitor->after_handle_event(m_currentEvent,
					                              m_status_line_insert_hint,
					                              m_is_event_on_above))
      {
        // It is possible to deallocate the event:
        deallocate_event(m_currentEvent);
      }

      // We are done with the current event - remove it from the queue.
      m_queue->erase(eventIter);
      eventIter = m_queue->begin();
    }

    return;
  }


  /*! create Event object for each input point */
  template <class PointInputIterator>
  void _init_points(PointInputIterator points_begin,
                   PointInputIterator points_end,
                   Attribute type)
  {
    for(PointInputIterator iter = points_begin;
        iter != points_end;
        ++iter)
    {
      _init_point(*iter, type);
    }
  }


  /*! for each curve create a Subcurve object and two Event objects */
  template<class CurveInputIterator>
  void _init_curves(CurveInputIterator curves_begin,
                    CurveInputIterator curves_end)
  {
    unsigned int   index = 0;
    for(CurveInputIterator iter = curves_begin;
        iter != curves_end;
        ++iter, ++index)
    {
      _init_curve(*iter, index);
    }
  }


  /*! Init the sweep algorithm */
  template<class CurveInputIterator>
  void _init_sweep(CurveInputIterator curves_begin,
                  CurveInputIterator curves_end)
  {
    m_num_of_subCurves = std::distance(curves_begin, curves_end);

    _init_structures();

    //init the curves
    _init_curves(curves_begin, curves_end);
  }


  /*! Init the data structures for the sweep algorithm */
  virtual void _init_structures()
  {
    CGAL_assertion(m_queue->empty() && (m_statusLine.size() == 0));

    
     //allocate all of the Subcure objects as one block
    m_subCurves = m_subCurveAlloc.allocate(m_num_of_subCurves);

  }


  /*! Compete the sweep (compete data strcures) */
  virtual void _complete_sweep()
  {
    CGAL_assertion(m_queue->empty() && (m_statusLine.size() == 0));

     for(unsigned int i=0 ; i < m_num_of_subCurves; ++i)
    m_subCurveAlloc.destroy(m_subCurves+i);

    if(m_num_of_subCurves) //if its zero, nothing to deallocate
      m_subCurveAlloc.deallocate(m_subCurves,m_num_of_subCurves);
  }


  /*!
   * Initialize an event associated with a point.
   * \param p The given point.
   */
  void _init_point(const Point_2& pt, Attribute type)
  {
    const std::pair<Event*, bool>& pair_res = push_event(pt, type);
    if(! pair_res.second)
      m_visitor-> update_event(pair_res.first, pt);

    m_visitor -> init_event(pair_res.first);
  }
  

  /*!
   * Initialize an event associated with an x-monotone curve.
   * \param curve The given x-monotone curve.
   * \param indec Its unique index.
   */
  void _init_curve(const X_monotone_curve_2 &curve,unsigned int index)
  {
     // construct a Subcurve object
    m_subCurveAlloc.construct(m_subCurves+index, m_masterSubcurve);

    (m_subCurves+index)->init(curve);

    _init_endpoint(curve, MAX_END, m_subCurves+index);
    _init_endpoint(curve, MIN_END, m_subCurves+index);
     
    return;
  }

  Event* _init_endpoint (const X_monotone_curve_2& cv,
                         Curve_end ind,
                         Subcurve* sc)
  {
    Attribute end_attr = 
      (ind == MIN_END) ? Base_event::LEFT_END : Base_event::RIGHT_END;
    Point_2 end_point;

    Boundary_type inf_x = m_traits->boundary_in_x_2_object()(cv, ind);
    Boundary_type inf_y = m_traits->boundary_in_y_2_object()(cv, ind);
    
    std::pair<Event*, bool> pair_res;
    if(inf_x != NO_BOUNDARY)
    {
      pair_res =
        push_event(cv, end_attr , inf_x, inf_y, ind, sc);
    }
    else
    {
      if(inf_y != NO_BOUNDARY)
      {
        pair_res =
          push_event(cv, end_attr , inf_x, inf_y, ind, sc);
      }
      else
      {
        end_point = 
          ((ind == MIN_END) ? m_traits->construct_min_vertex_2_object()(cv):
                              m_traits->construct_max_vertex_2_object()(cv));
                              
        pair_res =
          push_event(end_point, end_attr, sc);

      }
    }
    
    Event* e = pair_res.first;
    if(pair_res.second == false) // the event already exist
    {
      m_visitor ->update_event(e, end_point, cv , false); //new notification function!!
    }

    return e;
  }
  

  /*! Handle the subcurve to the left of the current event point. */
  virtual void _handle_left_curves()
  { 
    CGAL_PRINT("Handling left curve" << std::endl;);

    m_is_event_on_above = false;

    if(! m_currentEvent->has_left_curves())
    { 
      _handle_event_without_left_curves();
      if(m_currentEvent->is_finite())
      {
        if(m_is_event_on_above)
        {
          // current event is on the interior of existing curve at the Y-str,
          // it can allowed only if the event is an isolated query point
          CGAL_assertion(!m_currentEvent -> has_right_curves() &&
                          m_currentEvent -> is_query());

          //m_is_event_on_above = true;
          m_visitor->before_handle_event(m_currentEvent);
        }
        else
          m_visitor->before_handle_event(m_currentEvent);
      }
      else
        m_visitor->before_handle_event(m_currentEvent);
       
      //nothing else to do (no left curves)
      return;
    }
    
        

    CGAL_PRINT("left curves before sorting: "<<"\n";);
    CGAL_SL_DEBUG(if (m_currentEvent->left_curves_begin() != 
		 m_currentEvent->left_curves_end() )
             {
               m_currentEvent->Print();
             });
    // determine the order of left curves by the Y-structure
    _sort_left_curves();
    m_visitor->before_handle_event(m_currentEvent);

    CGAL_PRINT("left curves after sorting: "<<"\n";);
    CGAL_SL_DEBUG(if (m_currentEvent->left_curves_begin() != 
		 m_currentEvent->left_curves_end() )
             {
               m_currentEvent->Print();
             });

    EventCurveIter left_iter = m_currentEvent->left_curves_begin();
    while(left_iter != m_currentEvent->left_curves_end())
    {
      Subcurve *leftCurve = *left_iter; 
      m_visitor->add_subcurve(leftCurve->get_last_curve(), leftCurve);
      ++left_iter;

      //remove curve from the status line 
      _remove_curve_from_status_line(leftCurve);    
    }
    CGAL_PRINT( "Handling left curve END" << std::endl;);
      
    return;
  }

  void _handle_event_without_left_curves()
  {
     if(m_currentEvent->is_finite())
     {
       const std::pair<StatusLineIter, bool>& pair_res =
         m_statusLine.find_lower (m_currentEvent->point(), 
				                          m_statusLineCurveLess);

       m_status_line_insert_hint = pair_res.first;
       m_is_event_on_above = pair_res.second;

       return;
     }
        
    // its an event at infinity
    Boundary_type inf_x = m_currentEvent->infinity_at_x();
    if(inf_x == MINUS_INFINITY)
      m_status_line_insert_hint = m_statusLine.end();
    else
    {
      CGAL_assertion(inf_x != PLUS_INFINITY); //event at plus infinity x
                                              // must have left curve
      Boundary_type inf_y = m_currentEvent->infinity_at_y();
      if(inf_y == MINUS_INFINITY)
        m_status_line_insert_hint = m_statusLine.begin();
      else
      {
        CGAL_assertion(inf_y == PLUS_INFINITY);
        m_status_line_insert_hint = m_statusLine.end();
      }
    }
      
    return;
  }
  /*!
   * Sort the left subcurves of an event point according to their order in
   * their status line (no geometric comprasions are needed).
   */
  void _sort_left_curves()
  {
    CGAL_assertion(m_currentEvent->has_left_curves());
   
    Subcurve *curve = *(m_currentEvent->left_curves_begin());
    StatusLineIter slIter = curve->get_hint();
    CGAL_assertion(*slIter == curve);
   
    
    for (++slIter; slIter != m_statusLine.end(); ++slIter)
    {
      if(std::find(m_currentEvent->left_curves_begin(),
                   m_currentEvent->left_curves_end(),
                   *slIter) ==
                   m_currentEvent->left_curves_end())
         break;
    }
    StatusLineIter end (slIter);

    slIter = curve->get_hint();
    if(slIter == m_statusLine.begin())
    {
      m_currentEvent->replace_left_curves(slIter,end);
      return;
    }
    --slIter;
    for(;slIter != m_statusLine.begin(); --slIter)
    {
      if( std::find(m_currentEvent->left_curves_begin(),
                    m_currentEvent->left_curves_end(),
                    *slIter) == m_currentEvent->left_curves_end())
      {
        m_currentEvent->replace_left_curves(++slIter,end);
        return;
      }
    }
    if(std::find(m_currentEvent->left_curves_begin(),
                  m_currentEvent->left_curves_end(),
                  *slIter) == m_currentEvent->left_curves_end())
    {
      m_currentEvent->replace_left_curves(++slIter,end);;
    }
    else
    {
      m_currentEvent->replace_left_curves(slIter,end);
    }
  }

 


  /*! Handle the subcurve to the left of the current event point. */
  virtual void _handle_right_curves()
  {
    CGAL_PRINT("Handling right curves (" ;)
    CGAL_PRINT(m_currentEvent->point() << ")\n";)
    
    if(! m_currentEvent->has_right_curves())
      return;


    // Loop over the curves to the right of the current event and handle them:
    // since we are at allways the beginning of a curve, we just insert 
    // it to the status line.
   
    EventCurveIter currentOne = m_currentEvent->right_curves_begin();
    EventCurveIter rightCurveEnd = m_currentEvent->right_curves_end();

    
    while ( currentOne != rightCurveEnd )
    {
      CGAL_PRINT_INSERT(*currentOne);
      StatusLineIter slIter = 
        m_statusLine.insert_before(m_status_line_insert_hint, *currentOne);
      ((Subcurve*)(*currentOne))->set_hint(slIter);
        
      CGAL_SL_DEBUG(PrintStatusLine(););
      ++currentOne;
    }        
      
    CGAL_SL_DEBUG(PrintStatusLine(););
  }


  /*!
   * Add a subcurve to the right of an event point.
   * \param event The event point.
   * \param curve The subcurve to add.
   * \return (true) if an overlap occured; (false) otherwise.
   */
  virtual bool _add_curve_to_right (Event* event, Subcurve* curve,
                                    bool /* overlap_exist*/ = false)
  {
    std::pair<bool, EventCurveIter> pair_res = 
      event->add_curve_to_right(curve, m_traits);

    CGAL_assertion(!pair_res.first);
      return (false);
  }

  


  /*! Remove a curve from the status line. */
  void _remove_curve_from_status_line (Subcurve *leftCurve);

#ifdef CGAL_DEBUG_SWEEP_LINE
  void PrintEventQueue();
  void PrintSubCurves();
  void PrintStatusLine();
  void PrintInfinityType(Boundary_type x, Boundary_type y);
  void PrintEvent(const Event* e);
#endif

protected:

  /*! a  traits object */
  Traits_adaptor *m_traits;

  /*! indicates if the traits object was allocated by the sweep */
  bool m_traitsOwner;

   /*! a pointer to the current event */
  Event *m_currentEvent;

  /*! Y-str comprasion functor */
  StatusLineCurveLess m_statusLineCurveLess;

  EventLess   m_queueEventLess;
  /*! The queue of events (intersection points) to handle */
  EventQueue *m_queue;

  /*! The events that have been allocated (an not freed). */
  Allocated_events_set    m_allocated_events;

  /*! The subcurves array */
  Subcurve *m_subCurves;

  /*! The status line (Y-str) */
  StatusLine m_statusLine;
 
  /*! An iterator of the  status line that is used as a hint for inserts. */
  StatusLineIter m_status_line_insert_hint;

  /*! indicates if current event is on the interior of existing curve, it may
   *  happen only with events that are associated with isolated query points
   */
  bool m_is_event_on_above;

  /*! An allocator for the events objects */
  EventAlloc m_eventAlloc;

  /*! An allocator for the Subcurve objects */
  SubCurveAlloc m_subCurveAlloc;

  /*! a master Event (created once by the constructor) for the allocator's 
   *  usgae. */
  Event m_masterEvent;

  /*! a master Subcurve (created once by the constructor) for the allocator's 
   *  usgae. */
  Subcurve m_masterSubcurve;

  /*! The num of subcurves  */
  unsigned int m_num_of_subCurves;

  /*! a pointer to the visitor object which will be notidifed during sweep */
  SweepVisitor* m_visitor;

 

  /*! Allocate an event object */
  Event* _allocate_event(const Point_2& pt, Attribute type,Boundary_type,Boundary_type)
  {
    Event *e =  m_eventAlloc.allocate(1); 
    m_eventAlloc.construct(e, m_masterEvent);
    e->init(pt, type);

    m_allocated_events.insert(e);

    return e;
  }

    Event* _allocate_event(const Point_2& pt, Attribute type)
    {
       Event *e =  m_eventAlloc.allocate(1);
       m_eventAlloc.construct(e, m_masterEvent);
       e->init(pt, type);

       m_allocated_events.insert(e);
       return e;
    }

  /*! Push a finite event point to x-structure (m_queue) iff it doesnt exist */
  std::pair<Event*, bool> push_event(const Point_2& pt,
                                     Attribute type,
                                     Subcurve* sc = NULL)
  {
    Event*    e;  
    
    const std::pair<EventQueueIter, bool>& pair_res =
      m_queue->find_lower(pt, m_queueEventLess);
    bool exist = pair_res.second;

    if (! exist)
    {
      // We have a new event
      e = _allocate_event(pt, type);
      e->set_finite();

      if(sc != NULL)
      {
        if(type == Base_event::LEFT_END)
        {
          sc->set_left_event(e);
          _add_curve_to_right(e, sc);
        }
        else
        {
          CGAL_assertion(type == Base_event::RIGHT_END);
          sc->set_right_event(e);
          e->add_curve_to_left(sc);
        }
      }

      m_queue->insert_before(pair_res.first, e);
    }
    else
    {
      // The event already exsits
      e = *(pair_res.first);
      e->set_attribute(type);
      e->set_finite();
      if(sc != NULL)
      {
        if(type == Base_event::LEFT_END)
        {
          sc->set_left_event(e);
          _add_curve_to_right(e, sc);
        }
        else
        {
          CGAL_assertion(type == Base_event::RIGHT_END);
          sc->set_right_event(e);
          e->add_curve_to_left(sc);
        }
      }
    }
    CGAL_PRINT_NEW_EVENT(pt, e);
    
    return (std::make_pair(e, !exist));
  }

  /*! Push an event at infinity to x-structure (m_queue) iff it doesnt exist (overlap) */
  std::pair<Event*, bool> push_event(const X_monotone_curve_2& cv,
                                     Attribute type,
                                     Boundary_type x_inf,
                                     Boundary_type y_inf = NO_BOUNDARY,
                                     Curve_end     ind = MIN_END,
                                     Subcurve* sc = NULL)
  {
    Event*    e;  
    
    m_queueEventLess.set_boundary_in_x(x_inf);
    m_queueEventLess.set_boundary_in_y(y_inf);
    m_queueEventLess.set_index(ind);

    const std::pair<EventQueueIter, bool>& pair_res =
      m_queue->find_lower(cv, m_queueEventLess);
    bool exist = pair_res.second;

    if (! exist)
    {
      // We have a new event
      Point_2 pt = Point_2();
      e = _allocate_event(pt, type);
      _set_attributes_of_infinity(e, x_inf, y_inf);
      if(sc != NULL)
      {
        if(type == Base_event::LEFT_END)
        {
          sc->set_left_event(e);
          _add_curve_to_right(e, sc);
        }
        else
        {
          sc->set_right_event(e);
          e->add_curve_to_left(sc);
        }
      }
      m_queue->insert_before(pair_res.first, e);
    }
    else
    {
      // The event already exsits
      e = *(pair_res.first);
      e->set_attribute(type);
      _set_attributes_of_infinity(e, x_inf, y_inf);
      if(sc != NULL)
      {
        if(type == Base_event::LEFT_END)
        {
          sc->set_left_event(e);
          _add_curve_to_right(e, sc);
        }
        else
        {
          sc->set_right_event(e);
          e->add_curve_to_left(sc);
        }
      }
    }
    
    return (std::make_pair(e, !exist));
  }

  void _set_attributes_of_infinity(Event* e,
                                   Boundary_type x_inf,
                                   Boundary_type y_inf)
  {
    if(x_inf == MINUS_INFINITY)
      e->set_minus_infinite_x();
    else
      if(x_inf == PLUS_INFINITY)
        e->set_plus_infinite_x();
      else
        e->set_finite_x();
        
    if(y_inf == MINUS_INFINITY)
          e->set_minus_infinite_y();
    else
      if(y_inf == PLUS_INFINITY)
        e->set_plus_infinite_y();
      else
        e->set_finite_y();
  }

};



/*!
 * Remove a curve from the status line for good.
 *
 * @param leftCurve a pointer to the curve that is about to be deleted
 * @return 
 */
template <class Traits_,
          class SweepVisitor,
          class CurveWrap,
          class SweepEvent,
          typename Allocator>
inline void Basic_sweep_line_2<Traits_,
                               SweepVisitor,                               
                               CurveWrap,
                               SweepEvent,
                               Allocator>::
_remove_curve_from_status_line(Subcurve *leftCurve)
                              
{
  CGAL_PRINT("remove_curve_from_status_line\n";);
  CGAL_SL_DEBUG(PrintStatusLine(););
  CGAL_SL_DEBUG(leftCurve->Print(););

  StatusLineIter sliter = leftCurve->get_hint(); 
  m_status_line_insert_hint = sliter; ++m_status_line_insert_hint; 

  CGAL_assertion(sliter!=m_statusLine.end());
  m_statusLine.erase(sliter);
  CGAL_PRINT("remove_curve_from_status_line Done\n";)
} 



//DEBUG UTILITIES
#ifdef CGAL_DEBUG_SWEEP_LINE
  #include <CGAL/Sweep_line_2/Sweep_line_2_debug.h>
#endif

CGAL_END_NAMESPACE

#endif
