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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 (based on old version by Tali Zvi)

#ifndef CGAL_BASIC_SWEEP_LINE_2_H
#define CGAL_BASIC_SWEEP_LINE_2_H

#include <map>
#include <vector>
#include <algorithm>
#include <iterator>
#include <CGAL/assertions.h>
#include <CGAL/memory.h>
#include <CGAL/Sweep_line_2/Sweep_line_functors.h>
#include <CGAL/Sweep_line_2/Sweep_line_subcurve.h>
#include <CGAL/Sweep_line_2/Sweep_line_event.h>
#include <CGAL/Multiset.h>
//#include <CGAL/Sweep_line_2/Sweep_line_rb_tree.h>


#ifndef VERBOSE

#define SL_DEBUG(a)
#define PRINT_INSERT(a)
#define PRINT_ERASE(a)
#define PRINT_NEW_EVENT(p, e) 
#define PRINT(a)


#else

#define SL_DEBUG(a) {a}
#define PRINT_INSERT(a) { std::cout << "+++ inserting "; \
                          (a)->Print(); \
                          std::cout << "    currentPos = "  \
                                    << m_currentEvent->get_point() \
                                    << "\n"; \
                          }
#define PRINT_ERASE(a)  { std::cout << "--- erasing " ; \
                          (a)->Print(); }
#define PRINT_NEW_EVENT(p, e) \
{ std::cout << "%%% a new event was created at " << (p) << std::endl; \
  (e)->Print(); }
#define PRINT(a) { std::cout << a ; }

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
  typedef typename Traits::Point_2                       Point_2;
  typedef typename Traits::X_monotone_curve_2            X_monotone_curve_2;

  typedef SweepEvent                                     Event;
  typedef Point_less_functor<Traits>                     PointLess;
  typedef std::map<Point_2 , Event*, PointLess>          EventQueue; 
  typedef typename EventQueue::iterator                  EventQueueIter;
  typedef typename EventQueue::value_type                EventQueueValueType;

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

 

  /*!
   * Constructor.
   * \param visitor A pointer to a sweep-line visitor object.
   */
  Basic_sweep_line_2 (SweepVisitor* visitor) :
      m_traits(new Traits()),
      m_traitsOwner(true),
      m_statusLineCurveLess(m_traits, &m_currentEvent),
      m_queue(new EventQueue(PointLess(m_traits))),
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
      m_traits(traits),
      m_traitsOwner(false),
      m_statusLineCurveLess(m_traits, &m_currentEvent),
      m_queue(new EventQueue(PointLess(m_traits))),
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
      this ->deallocate_event(qiter->second);
    }
    this -> m_statusLine.clear();
    m_status_line_insert_hint = this -> m_statusLine.begin();
  
    CGAL_assertion(!m_queue->empty());
    EventQueueIter second = m_queue->begin(); ++second;
    m_queue->erase(second, m_queue->end());
  }


   /*! Deallocate event object, it is a public method to allow the visitor
    *  to manage the events deallocation (if he wants to) 
    */
  void deallocate_event(Event* event)
  {
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
      m_currentEvent = eventIter->second;

      PRINT("------------- " 
            << m_currentEvent->get_point() 
            << " --------------"
            << std::endl;);
      SL_DEBUG(PrintStatusLine();
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
    const Point_2 &left_end =
      m_traits->construct_min_vertex_2_object()(curve);
    const Point_2 &right_end =
      m_traits->construct_max_vertex_2_object()(curve);

    
    // Handle the right endpoint of the curve.
    const std::pair<Event*, bool>& pair_res1 =
      push_event(right_end, Base_event::RIGHT_END);

    Event* right_event = pair_res1.first;
    if(pair_res1.second == false) // the event already exist
    {
      m_visitor ->update_event(right_event, right_end, curve , false); //new notification function!!
    }
    
    // Handle the left endpoint of the curve.
    const std::pair<Event*, bool>& pair_res2 =
      push_event(left_end, Base_event::LEFT_END); 

    Event* left_event = pair_res2.first;
    if(pair_res2.second == false) // the even already exist
    {
      m_visitor ->update_event(left_event, left_end, curve, true); //new notification function!!
    }
    
    // construct a Subcurve object
    m_subCurveAlloc.construct(m_subCurves+index, m_masterSubcurve);

    (m_subCurves+index)->init(curve, left_event, right_event);
    
    right_event->add_curve_to_left(m_subCurves+index);  
    _add_curve_to_right(left_event, m_subCurves+index);
    
    return;
  }
  

  /*! Handle the subcurve to the left of the current event point. */
  virtual void _handle_left_curves()
  { 
    PRINT("Handling left curve" << std::endl;);

    m_is_event_on_above = false;

    if(! m_currentEvent->has_left_curves())
    { 
      m_statusLineCurveLess.set_is_equal(false);
      m_status_line_insert_hint =
        m_statusLine.lower_bound (m_currentEvent->get_point(), 
				                          m_statusLineCurveLess);

      m_is_event_on_above = m_statusLineCurveLess.is_equal();

      if(m_is_event_on_above)
      {
        // current event is on the interior of existing curve at the Y-str,
        // it can allowed only if the event is an isolated query point
        CGAL_assertion(!m_currentEvent -> has_right_curves() &&
                        m_currentEvent -> is_query());

         m_is_event_on_above = true;
         m_visitor->before_handle_event(m_currentEvent);
      }
      else
        m_visitor->before_handle_event(m_currentEvent);
       
      //nothing else to do (no left curves)
      return;
    }
    
        

    PRINT("left curves before sorting: "<<"\n";);
    SL_DEBUG(if (m_currentEvent->left_curves_begin() != 
		 m_currentEvent->left_curves_end() )
             {
               m_currentEvent->Print();
             });
    // determine the order of left curves by the Y-structure
    _sort_left_curves();
    m_visitor->before_handle_event(m_currentEvent);

    PRINT("left curves after sorting: "<<"\n";);
    SL_DEBUG(if (m_currentEvent->left_curves_begin() != 
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
    PRINT( "Handling left curve END" << std::endl;);
      
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
    PRINT("Handling right curves (" ;)
    PRINT(m_currentEvent->get_point() << ")\n";)
    
    if(! m_currentEvent->has_right_curves())
      return;


    // Loop over the curves to the right of the current event and handle them:
    // since we are at allways the beginning of a curve, we just insert 
    // it to the status line.
   
    EventCurveIter currentOne = m_currentEvent->right_curves_begin();
    EventCurveIter rightCurveEnd = m_currentEvent->right_curves_end();

    
    while ( currentOne != rightCurveEnd )
    {
      PRINT_INSERT(*currentOne);
      StatusLineIter slIter = 
        m_statusLine.insert_before(m_status_line_insert_hint, *currentOne);
      ((Subcurve*)(*currentOne))->set_hint(slIter);
        
      SL_DEBUG(PrintStatusLine(););
      ++currentOne;
    }        
      
    SL_DEBUG(PrintStatusLine(););
  }


  /*!
   * Add a subcurve to the right of an event point.
   * \param event The event point.
   * \param curve The subcurve to add.
   * \return (true) if an overlap occured; (false) otherwise.
   */
  virtual bool _add_curve_to_right (Event* event, Subcurve* curve,
                                    bool overlap_exist = false)
  {
    std::pair<bool, EventCurveIter> pair_res = 
      event->add_curve_to_right(curve, m_traits);

    CGAL_assertion(!pair_res.first);
      return (false);
  }

  


  /*! Remove a curve from the status line. */
  void _remove_curve_from_status_line (Subcurve *leftCurve);

#ifdef VERBOSE
  void PrintEventQueue();
  void PrintSubCurves();
  void PrintStatusLine();
#endif

protected:

  /*! a  traits object */
  Traits *m_traits;

  /*! indicates if the traits object was allocated by the sweep */
  bool m_traitsOwner;

   /*! a pointer to the current event */
  Event *m_currentEvent;

  /*! Y-str comprasion functor */
  StatusLineCurveLess m_statusLineCurveLess;

  /*! the queue of events (intersection points) to handle */
  EventQueue *m_queue;

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
  Event* allocate_event(const Point_2& pt, Attribute type)
  {
    Event *e =  m_eventAlloc.allocate(1); 
    m_eventAlloc.construct(e, m_masterEvent);
    e->init(pt, type);
    return e;
  }

  /*! Push event point to x-structure (m_queue) iff it doesnt exist */
  std::pair<Event*, bool> push_event(const Point_2& pt, Attribute type)
  {
    Event*    e;  
    const std::pair<EventQueueIter, bool>& insertion_res =
      m_queue->insert(EventQueueValueType(pt,0));
    bool inserted = insertion_res.second;

    if (inserted == true)
    {
      // We have a new event
      e = allocate_event(pt, type);
      (insertion_res.first)->second = e;
    }
    else
    {
      // The event already exsits
      e = (insertion_res.first)->second;
      e->set_attribute(type);
    }
    PRINT_NEW_EVENT(pt, e);
    return (std::make_pair(e, inserted));
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
  PRINT("remove_curve_from_status_line\n";);
  SL_DEBUG(PrintStatusLine(););
  SL_DEBUG(leftCurve->Print(););

  StatusLineIter sliter = leftCurve->get_hint(); 
  m_status_line_insert_hint = sliter; ++m_status_line_insert_hint; 

  CGAL_assertion(sliter!=m_statusLine.end());
  m_statusLine.erase(sliter);
  PRINT("remove_curve_from_status_line Done\n";)
} 



//DEBUG UTILITIES
#ifdef VERBOSE
  #include <CGAL/Sweep_line_2/Sweep_line_2_debug.h>
#endif

CGAL_END_NAMESPACE

#endif
