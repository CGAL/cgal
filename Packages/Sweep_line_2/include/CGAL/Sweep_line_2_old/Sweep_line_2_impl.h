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
// Author(s)     : Tali Zvi <talizvi@post.tau.ac.il>,
//                 Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef CGAL_SWEEP_LINE_2_IMPL_H
#define CGAL_SWEEP_LINE_2_IMPL_H

#include <map>
#include <set>
#include <list>
#include <vector>
#include <algorithm>
#include <iterator>
#include <CGAL/assertions.h>
#include <CGAL/memory.h>
#include <CGAL/Object.h>
#include <CGAL/Sweep_line_2_old/Sweep_line_functors.h>
#include <CGAL/Sweep_line_2_old/Sweep_line_traits.h>
#include <CGAL/Sweep_line_2_old/Sweep_line_event.h>
#include <CGAL/Sweep_line_2_old/Sweep_line_subcurve.h>
#include <CGAL/Sweep_line_2_old/Sweep_line_rb_tree.h>



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
                          std::cout << "    currentPos = "  \
                                    << m_currentEvent->get_point() \
                                    << "\n"; \
                          }
#define PRINT_ERASE(a)  { std::cout << "--- erasing " ; \
                          (a)->Print(); }
#define PRINT_NEW_EVENT(p, e) \
{ std::cout << "%%% a new event was created at " << (p) << std::endl; \
  (e)->Print(); }
#define DBG(a) { std::cout << a << std::endl; }
#define STORE_RESULT(a) {a}
#endif



CGAL_BEGIN_NAMESPACE

/*!
  Sweep_line_2_impl is a class that implements the sweep line algorithm
  based on the algorithm of Bentley and Ottmann.
  It extends the algorithm to support not only segments but polylines and 
  general curves as well.
  The curves are defined by the traits class that is one of the template 
  arguments.

  The algorithm is also extended to support the following degenerate cases:
  - non x-monotone curves
  - vertical segments
  - multiple (more then two) segments intersecting at one point
  - curves beginning and ending on other curves.
  - overlapping curves

  An extension of this algorithm is used to produce a planar map containing
  the input curves efficiently. \sa Pmwx_aggregate_insert.

  General flow:
  After the initialization stage, the events are handled from left to right.

  For each event

    First pass - handles special cases in which curves start or end 
                 at the interior of another curve
    Handle left curves - iterate over the curves that intersect 
                 at the event point and defined to the left of the 
                 event. 
    Handle right curves - iterate over the curves that intersect 
                 the event point and defined to the right of the 
                 event point. This is where new intersection points 
                 are calculated.
  End

  Convensions through out the code:
  In order to make the code as readable as possible, some convensions were 
  made in regards to variable naming:

    xp - is the intersection point between two curves
    slIter - an iterator to the status line, always points to a curve.

*/
template < class SweepLineTraits_2,
           class SweepEvent,
           class CurveWrap,
           class SweepVisitor,
           typename Allocator = CGAL_ALLOCATOR(int) >
class Sweep_line_2_impl
{
public:

  typedef SweepLineTraits_2 Traits;
  typedef typename Traits::Point_2 Point_2;
  typedef typename Traits::Curve_2 Curve_2;
  typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;

  typedef SweepEvent Event;
  typedef Point_less_functor<Traits> PointLess;
  typedef std::map< Point_2 , Event*, PointLess> EventQueue; 
  typedef typename EventQueue::iterator EventQueueIter;
  typedef typename EventQueue::value_type EventQueueValueType;

  typedef typename Event::SubCurveIter EventCurveIter;

  typedef CurveWrap Subcurve;
  typedef std::list<Subcurve*> SubCurveList;
  typedef typename SubCurveList::iterator SubCurveListIter;
   
  typedef Status_line_curve_less_functor<Traits, Subcurve> StatusLineCurveLess;
  typedef Red_black_tree<Subcurve*, StatusLineCurveLess, CGAL_ALLOCATOR(int)> StatusLine;
  typedef typename StatusLine::iterator StatusLineIter;
 
  typedef typename Allocator::template rebind<Event>      EventAlloc_rebind;
  typedef typename EventAlloc_rebind::other               EventAlloc;

  typedef typename Allocator::template rebind<Subcurve>   SubcurveAlloc_rebind;
  typedef typename SubcurveAlloc_rebind::other            SubCurveAlloc;

 


  Sweep_line_2_impl(SweepVisitor* visitor) :
      //m_traits(new Traits()),
      m_sweep_line_traits(&m_traits),
      //m_traitsOwner(true),
      m_statusLineCurveLess(&m_traits),
      m_queue(new EventQueue(PointLess(&m_traits))),
      m_statusLine(m_statusLineCurveLess),
      m_xcurves(0),
      m_status_line_insert_hint(m_statusLine.begin()),
      m_num_of_subCurves(0),
#ifndef NDEBUG
      m_eventId(0),
#endif
      m_visitor(visitor)
  {
    m_visitor->attach(this);
  }


  Sweep_line_2_impl(Traits *t, SweepVisitor* visitor) :
      m_traits(*t),
      m_sweep_line_traits(&m_traits),
      //m_traitsOwner(false),
      m_statusLineCurveLess(&m_traits),
      m_queue(new EventQueue(PointLess(&m_traits))),
      m_statusLine(m_statusLineCurveLess),
      m_xcurves(0),
      m_status_line_insert_hint(m_statusLine.begin()),
      m_num_of_subCurves(0),
#ifndef NDEBUG
      m_eventId(0),
#endif
      m_visitor(visitor)
  {
    m_visitor->attach(this);
  }


  virtual ~Sweep_line_2_impl();

 


  public:


  /*
   *  - x-monotonize the input curves
   *  - for each end point of each curve create an event
   *  - initialize the event queue
   *  -
   */
  template<class CurveInputIterator>
  void init(CurveInputIterator curves_begin, CurveInputIterator curves_end)
  {
    //std::cout<<"distance is: " <<std::distance(curves_begin,curves_end) <<"!!\n";
    //m_xcurves.reserve(std::distance(curves_begin,curves_end));
    for(CurveInputIterator iter = curves_begin; iter != curves_end; ++iter)
    {
      m_traits.curve_make_x_monotone(*iter, std::back_inserter(m_xcurves));
    }
    m_num_of_subCurves = m_xcurves.size();
    m_subCurves = m_subCurveAlloc.allocate(m_num_of_subCurves);

    for(unsigned int index = 0; index < m_num_of_subCurves ; ++index)
    {
      init_curve(m_xcurves[index],index);
    }
  }

  template<class CurveInputIterator, class XCurveInputIterator>
  void init(CurveInputIterator curves_begin, CurveInputIterator curves_end,
            XCurveInputIterator xcurves_begin, XCurveInputIterator xcurves_end)
  {
   /* std::vector<X_monotone_curve_2> m_xcurves;
    m_xcurves.reserve(std::distance(curves_begin,curves_end) +
                      std::distance(xcurves_begin,xcurves_end));  */                                                           
    for(CurveInputIterator iter = curves_begin; iter != curves_end; ++iter)
    {
      m_traits.curve_make_x_monotone(*iter, std::back_inserter(m_xcurves));
    }
    std::copy(xcurves_begin, xcurves_end, std::back_inserter(m_xcurves));
    m_num_of_subCurves = m_xcurves.size();
    m_subCurves = m_subCurveAlloc.allocate(m_num_of_subCurves);
    for(unsigned int index = 0; index < m_num_of_subCurves ; ++index)
    {
      init_curve(m_xcurves[index],index);
    }
  }

  template<class CurveInputIterator, class PointInputIterator>
  void init(CurveInputIterator curves_begin, CurveInputIterator curves_end,
            PointInputIterator points_begin, PointInputIterator points_end,bool)
  {
    std::copy(curves_begin, curves_end, std::back_inserter(m_xcurves));
    m_num_of_subCurves = m_xcurves.size();
    m_subCurves = m_subCurveAlloc.allocate(m_num_of_subCurves);
    for(unsigned int index = 0; index < m_num_of_subCurves ; ++index)
    {
      init_curve(m_xcurves[index],index);
    }
    for(PointInputIterator itr = points_begin; itr != points_end; ++itr)
    {
      init_point(*itr);
    }
  }




  void init_point(const Point_2& pt)
  {
    Event *e;
    //handle the right end of the curve
    const std::pair<EventQueueIter, bool>& insertion_res =
      m_queue->insert(EventQueueValueType(pt,0));  //the insertion return value

    if(insertion_res.second == true)  // its a new event
    {
      e = allocate_event(pt);
      m_visitor -> init_event(e);
     
    #ifndef NDEBUG
      e->id = m_eventId++;
    #endif
     
      (insertion_res.first)->second = e;
    }

    else // the event is already exist
    {
      e = (insertion_res.first)->second;
      m_visitor -> init_event(e);
    }
  }

  
  void init_curve(X_monotone_curve_2 &curve,unsigned int index);

  /*! The main loop to calculate intersections among the curves
   *  Looping over the events in the queue, for each event we first
   *  handle the curves that are to the left of the event point (i.e., 
   *  curves that we are done with), and then we look at the curves 
   *  to the right of the point, which means we attempt to find intersections
   *  between them and their neighbours on the sweep line.
   */
  void sweep()
  {
    EventQueueIter eventIter = m_queue->begin();
    while ( eventIter != m_queue->end() )
    {
      m_currentEvent = eventIter->second;

      SL_DEBUG(std::cout << "------------- " 
                         << m_currentEvent->get_point() 
                         << " --------------"
                         << std::endl;
               PrintStatusLine();
               m_currentEvent->Print(););
      
      handle_left_curves();
      m_queue->erase(eventIter);

      handle_right_curves();
      if(m_visitor->after_handle_event(m_currentEvent,
                                       m_status_line_insert_hint,
                                       m_is_event_on_above))
        deallocate_event(m_currentEvent);
      eventIter = m_queue->begin();
    }
  }

  StatusLineIter  StatusLine_begin()
  {
    return m_statusLine.begin();
  }

  StatusLineIter  StatusLine_end()
  {
    return m_statusLine.end();
  }


  protected:
  
  /*! For each left-curve, if it is the "last" subcurve, i.e., the 
   * event point is the right-edge of the original curve, the 
   * last sub curve is created and added to the result. Otherwise
   * the curve is added as is to the result.
   */
  void handle_left_curves()
  {
    //m_visitor->before_handle_event(m_currentEvent);
    
    SL_DEBUG(std::cout << "Handling left curve" << std::endl;);
    SL_DEBUG(if( m_currentEvent->left_curves_begin() != m_currentEvent->left_curves_end() )
             {
               m_currentEvent->Print();
             });
    if(!m_currentEvent -> has_left_curves())
    {
      /* this block takes care of
      //
      //           /
      //          /
      //       --------
      //          \
      //           \
      */
      SL_DEBUG(std::cout << " - handling special case " << std::endl;)
                                                     
      const std::pair<StatusLineIter, bool>& res =
        m_statusLine.lower_bound(m_currentEvent->get_point(), m_statusLineCurveLess);
      m_status_line_insert_hint = res.first;

      if(res.second) // indicates that current event point starts at the 
                     //interior of a curve at the y-str 
                     //(can also indicates overlap)
      {
        if(!m_currentEvent -> has_right_curves())
        {
          m_is_event_on_above = true;
          m_visitor->before_handle_event(m_currentEvent);
          return;
        }
        m_is_event_on_above = false;
        Subcurve* sc = *m_status_line_insert_hint;
        X_monotone_curve_2 last_curve = sc->get_last_curve();
        m_currentEvent->add_curve_to_left(sc); 
        bool is_overlap = add_curve_to_right(m_currentEvent, sc);
        X_monotone_curve_2 a,b;
        if ( sc->is_source_left_to_target() )
        {
          m_traits.curve_split(last_curve, a, b,m_currentEvent->get_point());
        }
        else
        {
          m_traits.curve_split(last_curve, b, a, m_currentEvent->get_point());
        }
        ++m_status_line_insert_hint; 
        m_statusLine.remove_at(res.first);
        if(!is_overlap)
        {
          sc->set_last_curve(b);
        }
        m_visitor->before_handle_event(m_currentEvent);
        m_visitor->add_subcurve(a, sc);
        return;
      }
      m_is_event_on_above = false;
      m_visitor->before_handle_event(m_currentEvent);
      return;
    }
        
    m_is_event_on_above = false; 
    sort_left_curves();
    m_visitor->before_handle_event(m_currentEvent);

     // indicates if the curve will be removed for good
    bool remove_for_good = false; 

    EventCurveIter left_iter = m_currentEvent->left_curves_begin();
    while(left_iter != m_currentEvent->left_curves_end())
    {
      Subcurve *leftCurve = *left_iter; 
    
      if((Event*)leftCurve->get_right_event() == m_currentEvent)
      {  
        remove_for_good = true;
        m_visitor->add_subcurve(leftCurve->get_last_curve(), leftCurve);

        if(leftCurve->get_orig_subcurve1() != NULL)
        {
          leftCurve->get_orig_subcurve1()->set_overlap_subcurve(NULL);
          leftCurve->get_orig_subcurve2()->set_overlap_subcurve(NULL);

          // clip the two subcurves according to end_overlap point
          Subcurve* res;
          if(( res = (Subcurve*) (leftCurve->get_orig_subcurve1() -> 
               clip(m_currentEvent))) != NULL)
          {
            add_curve_to_right(m_currentEvent, res);

          }

          if(( res = (Subcurve*)(leftCurve->get_orig_subcurve2() -> 
               clip(m_currentEvent))) != NULL)
          {
            add_curve_to_right(m_currentEvent, res);
          }
        }
      }
      else
      { 
        X_monotone_curve_2 a, b;
        const X_monotone_curve_2 &lastCurve = leftCurve->get_last_curve();
        if ( leftCurve->is_source_left_to_target() ) 
        {
          m_traits.curve_split(lastCurve, a, b, m_currentEvent->get_point());
          m_visitor->add_subcurve(a, leftCurve);
        }
        else
        {
          m_traits.curve_split(lastCurve, b, a, m_currentEvent->get_point());
          m_visitor->add_subcurve(a, leftCurve);
        }
        leftCurve->set_last_curve(b);
      }
       ++left_iter;

      // remove curve from the status line (also checks intersection 
      // between the neighbouring curves,only if the curve is removed for good)
      remove_curve_from_status_line(leftCurve, remove_for_good);    
    }
    SL_DEBUG(std::cout << "Handling left curve END" << std::endl;);
      
  }



 //sorts left curves according to the Y-str (no geometric comprasions)
  void sort_left_curves()
  {
    CGAL_assertion(m_currentEvent->has_left_curves());
    EventCurveIter leftCurveIter = m_currentEvent->left_curves_begin();
    while ( leftCurveIter != m_currentEvent->left_curves_end() )  
    {
      if( (*leftCurveIter)->get_overlap_subcurve() != NULL &&
           m_traits.compare_xy(
            (*leftCurveIter)->get_overlap_subcurve()->get_left_end(),
             m_currentEvent->get_point()) == SMALLER)
      {
        Subcurve *leftCurve = (Subcurve*)((*leftCurveIter)->getSubcurve());
        m_currentEvent->replace_right_curve((*leftCurveIter), leftCurve);
        *leftCurveIter = leftCurve;
      }
      ++leftCurveIter;
    }
    Subcurve *curve = *(m_currentEvent->left_curves_begin());
    StatusLineIter slIter = curve->get_hint();
    CGAL_assertion(*slIter == curve);
   
    
    for( ++slIter ;slIter != m_statusLine.end(); ++slIter)
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






  /*! Loop over the curves to the right of the status line and handle them:
   * - if we are at the beginning of the curve, we insert it to the status 
   *   line, then we look if it intersects any of its neighbours.
   * - if we are at an intersection point between two curves, we add them
   *   to the status line and attempt to intersect them with their neighbours. 
   * - We also check to see if the two intersect again to the right of the 
   *   point.
   */
  void handle_right_curves()
  {
    SL_DEBUG(std::cout << "Handling right curves (" ;)
    SL_DEBUG(std::cout << m_currentEvent->get_point() << ")\n";)
    
    if(! m_currentEvent->has_right_curves())
      return;

    int numRightCurves = m_currentEvent->get_num_right_curves();
    if(numRightCurves == 1)
    {
      SL_DEBUG(std::cout << " - beginning of curve " << std::endl;);
      SL_DEBUG(
        Subcurve *tmp1 = *(m_currentEvent->right_curves_begin());
        PRINT_INSERT(tmp1);
               );

      StatusLineIter slIter = 
        m_statusLine.insert_predecessor(m_status_line_insert_hint, 
                                *(m_currentEvent->right_curves_begin()));
      
      (*(m_currentEvent->right_curves_begin()))->set_hint(slIter); 
     
      SL_DEBUG(PrintStatusLine(););
     
      // if this is the only curve on the status line, nothing else to do
      if ( m_statusLine.size() == 1 )
        return;

      StatusLineIter prev = slIter;
      StatusLineIter next = slIter;
      ++next;
      if ( slIter != m_statusLine.begin() )
      {
        --prev;
        intersect(*prev, *(m_currentEvent->right_curves_begin()));
      }
      if ( next != m_statusLine.end() )
      { 
        intersect(*(m_currentEvent->right_curves_begin()), *next);
      } 
    }
    else  // numRightCurves > 1
    {
      EventCurveIter firstOne = m_currentEvent->right_curves_begin();
      EventCurveIter lastOne = m_currentEvent->right_curves_end(); --lastOne;
      EventCurveIter rightCurveEnd = m_currentEvent->right_curves_end();
      PRINT_INSERT(*firstOne);

      StatusLineIter slIter = m_statusLine.insert_predecessor(m_status_line_insert_hint, 
                                                    *firstOne);
      ((Subcurve*)(*firstOne))->set_hint(slIter);
        
      SL_DEBUG(PrintStatusLine(););
      if ( slIter != m_statusLine.begin() )
      { 
        //  get the previous curve in the y-str
        StatusLineIter prev = slIter; --prev;
        intersect(*prev, *slIter);
      }
      
      EventCurveIter currentOne = firstOne; ++currentOne;
      EventCurveIter prevOne = firstOne;
      while ( currentOne != rightCurveEnd )
      {
        PRINT_INSERT(*currentOne);
        slIter = m_statusLine.insert_predecessor(m_status_line_insert_hint, *currentOne);
        ((Subcurve*)(*currentOne))->set_hint(slIter);
          
        SL_DEBUG(PrintStatusLine(););
      
        intersect(*prevOne, *currentOne);
        prevOne = currentOne;
        ++currentOne;
      }        
      lastOne = currentOne; --lastOne;
        
      SL_DEBUG(PrintStatusLine(););
      StatusLineIter next = slIter; ++next;
      if ( next != m_statusLine.end() )
        intersect( *prevOne,*next);
    }
  }



  bool add_curve_to_right(Event* event, Subcurve* curve)
  {
    std::pair<bool, SubCurveListIter> pair_res = event->add_curve_to_right(curve);
    if(pair_res.first == true) //overlap
    {
      SL_DEBUG(std::cout<<"Overlap detected at right insertion...\n";);
      SubCurveListIter iter = pair_res.second;
      
       X_monotone_curve_2 overlap_cv;
       Object cv_obj =
        m_traits.nearest_intersection_to_right(curve->get_last_curve(),
                                                (*iter)->get_last_curve(), 
                                                event->get_point());
      if (CGAL::assign(overlap_cv, cv_obj))
      {
         //alocate new Subcure for the overlap
         Subcurve *overlap_sc = m_subCurveAlloc.allocate(1);
         m_subCurveAlloc.construct(overlap_sc,m_masterSubcurve);
         overlap_sc->init(overlap_cv);
         m_overlap_subCurves.push_back(overlap_sc);
         
         //get the right end of overlap_cv
         Point_2 end_overlap = m_traits.curve_target(overlap_cv);
         if(m_traits.compare_xy(end_overlap,
                                 m_traits.curve_source(overlap_cv)) == SMALLER)
         {
           end_overlap = m_traits.curve_source(overlap_cv);
         }
 
         //find the event assiciated with end_overlap point (right end point)
         EventQueueIter q_iter = m_queue->find( end_overlap );
         CGAL_assertion(q_iter!=m_queue->end());

         // set the members left event and right event of overlap_sc
         overlap_sc->set_left_event(event);
         overlap_sc->set_right_event((*q_iter).second);

         m_visitor->init_subcurve(overlap_sc);


         //remove curve, *iter from the left curves of end_overlap event
         ((*q_iter).second)->remove_curve_from_left(curve);
         ((*q_iter).second)->remove_curve_from_left(*iter);

         //add overlap_sc to the left curves
         ((*q_iter).second)->add_curve_to_left(overlap_sc);


         curve  -> set_overlap_subcurve(overlap_sc);
         (*iter)-> set_overlap_subcurve(overlap_sc);

         overlap_sc -> set_orig_subcurve1(*iter);
         overlap_sc -> set_orig_subcurve2(curve);  

          //replace current sub-curve (*iter) with the new sub-curve
         (*iter) = overlap_sc;
      }
    }
    return pair_res.first;
  }

  
  // utility methods 
  void intersect(Subcurve *c1, Subcurve *c2);
  void remove_curve_from_status_line(Subcurve *leftCurve,bool remove_for_good);

 
  
#ifndef NDEBUG
  void PrintEventQueue();
  void PrintSubCurves();
  void PrintStatusLine();
#endif

protected:

  /*! a pointer to a traits object */
  Traits m_traits;

  /*! an object that holds a static trait object 
   *  to be used by events and subcurves
   */
  Sweep_line_traits<Traits> m_sweep_line_traits;

  /*! Y-str comprasion functor */
  StatusLineCurveLess m_statusLineCurveLess;

  /*! the queue of events (intersection points) to handle */
  EventQueue *m_queue;

  /*! The subcurves array */
  Subcurve *m_subCurves;

  /*! The status line (Y-str) */
  StatusLine m_statusLine;

  /*! if non x-monotone are specified, this hold the x-monotone 
    curves created when splitting them into x-monotone curves. */
  std::vector<X_monotone_curve_2> m_xcurves;

  /*! a pointer to the current event */
  Event *m_currentEvent;

  /*! An iterator of the  status line that is used as a hint for inserts. */
  StatusLineIter m_status_line_insert_hint;

  /*! indicates if current event is on the interior of existing curve*/
  bool m_is_event_on_above;

  /*! An allocator for the events objects */
  EventAlloc m_eventAlloc;

  /*! An allocator for the Subcurve objects */
  SubCurveAlloc m_subCurveAlloc;

  /*! a master Event (created once by default c'tor)to be used by allocator */
  Event m_masterEvent;

  /*! a master Subcurve (created once by default c'tor) to be used by allocator */
  Subcurve m_masterSubcurve;

  /*! The num of subcurves  */
  unsigned int m_num_of_subCurves;

  /*! contains all of the new sub-curve creaed by overlap */
  SubCurveList m_overlap_subCurves;

#ifndef NDEBUG
  int m_eventId;
#endif

  /* pointer to the visitor object which will be notidifed during sweep */
  SweepVisitor* m_visitor;

 

  Event* allocate_event(const Point_2& pt)
  {
    Event *e =  m_eventAlloc.allocate(1); 
    m_eventAlloc.construct(e, m_masterEvent);
    e->init(pt);
    return e;
  }

  public:
  
  void deallocate_event(Event* event)
  {
    m_eventAlloc.destroy(event);
    m_eventAlloc.deallocate(event,1);
  }

};

template < class SweepLineTraits_2, class SweepEvent, class CurveWrap , class SweepVisitor,
          typename Allocator >
Sweep_line_2_impl<SweepLineTraits_2,SweepEvent,CurveWrap,SweepVisitor,
                  Allocator>::
~Sweep_line_2_impl() 
{  
  for(unsigned int i=0 ; i < m_num_of_subCurves; ++i)
    m_subCurveAlloc.destroy(m_subCurves+i);

  if(m_num_of_subCurves) //if its zero, nothing to deallocate
    m_subCurveAlloc.deallocate(m_subCurves,m_num_of_subCurves); // deallocate memory 

  for(SubCurveListIter itr = m_overlap_subCurves.begin();
      itr != m_overlap_subCurves.end();
      ++itr)
  {
    m_subCurveAlloc.destroy(*itr);
    m_subCurveAlloc.deallocate(*itr, 1);
  }
  delete m_queue;
}





/*! Given an x-monotone curve, create events for each end (if 
 *  one doesn't exist already). 
 *  For each curve create a Subcurve instance.
 */
template < class SweepLineTraits_2,
          class SweepEvent, class CurveWrap,class SweepVisitor,
          typename Allocator>
inline void 
Sweep_line_2_impl<SweepLineTraits_2,SweepEvent,CurveWrap,SweepVisitor,
                  Allocator>::
init_curve(X_monotone_curve_2 &curve,unsigned int j)
{ 
  Event *e ;
 
   m_subCurveAlloc.construct(m_subCurves+j,m_masterSubcurve);
  (m_subCurves+j)->init(curve);
 
  const Point_2 &left_end = ((m_subCurves+j)->is_source_left_to_target() ?
    m_traits.curve_source(curve) :
    m_traits.curve_target(curve));


  const Point_2 &right_end = ((m_subCurves+j)->is_source_left_to_target() ?
    m_traits.curve_target(curve) :
    m_traits.curve_source(curve));
    

    //handle the right end of the curve
    const std::pair<EventQueueIter, bool>& insertion_res =
      m_queue->insert(EventQueueValueType(right_end,0));  //the insertion return value

    if(insertion_res.second == true)  // its a new event
    {
      e = allocate_event(right_end);
     

    #ifndef NDEBUG
      e->id = m_eventId++;
    #endif
     
      (insertion_res.first)->second = e;
    }

    else // the event is already exist
    {
      e = (insertion_res.first)->second;
    }

    e->add_curve_to_left(m_subCurves+j);
    PRINT_NEW_EVENT(right_end, e);

    (m_subCurves+j)->set_right_event(e);


    // handle the left end of the curve
    const std::pair<EventQueueIter, bool>& insertion_res2 =
      m_queue->insert(EventQueueValueType(left_end,0));  //the insertion return value
     
    if(insertion_res2.second == true)
    {
      e = allocate_event(left_end);
    #ifndef NDEBUG
      e->id = m_eventId++;
    #endif
    
      (insertion_res2.first)->second = e;  
    }

    else // the event is already exist
    {
      e = (insertion_res2.first)->second;
    }
    add_curve_to_right(e, m_subCurves+j);
    PRINT_NEW_EVENT(left_end, e);

    (m_subCurves+j)->set_left_event(e);

    m_visitor->init_subcurve(m_subCurves+j);
  
}




/*!
 * When a curve is removed from the status line for good, its top and
 * bottom neighbors become neighbors. This method finds these cases and
 * looks for the intersection point, if one exists.
 *
 * @param leftCurve a pointer to the curve that is about to be deleted
 * @return an iterator to the position where the curve will be removed from.
 */
template <class SweepLineTraits_2,
          class SweepEvent, class CurveWrap, class SweepVisitor,
          typename Allocator>
inline void
Sweep_line_2_impl<SweepLineTraits_2,SweepEvent,CurveWrap,SweepVisitor,
                  Allocator>::
remove_curve_from_status_line(Subcurve *leftCurve, bool remove_for_good)
                              
{
  SL_DEBUG(std::cout << "remove_curve_from_status_line\n";);
  SL_DEBUG(PrintStatusLine(););
  SL_DEBUG(leftCurve->Print(););

  StatusLineIter sliter = leftCurve->get_hint(); 
  m_status_line_insert_hint = sliter; ++m_status_line_insert_hint;

  if(! remove_for_good)
  {
    m_statusLine.remove_at(sliter);
    SL_DEBUG(std::cout << "remove_curve_from_status_line Done\n";)
    return;
  }

  CGAL_assertion(sliter!=m_statusLine.end());
  if (sliter != m_statusLine.begin() && sliter != m_statusLine.maximum()) 
  {
    StatusLineIter prev = sliter; --prev;
    StatusLineIter next = sliter; ++next;
    
    // intersect *next with  *prev 
    intersect(*prev, *next);
  }
  m_statusLine.remove_at(sliter);
  SL_DEBUG(std::cout << "remove_curve_from_status_line Done\n";)
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


template < class SweepLineTraits_2,
           class SweepEvent, class CurveWrap,class SweepVisitor,
           typename Allocator >
inline void 
Sweep_line_2_impl<SweepLineTraits_2,SweepEvent,CurveWrap,SweepVisitor,
                  Allocator>::
intersect(Subcurve *c1, Subcurve *c2)
{
  SL_DEBUG(std::cout << "Looking for intersection between:\n\t";)
  SL_DEBUG(c1->Print();)
  SL_DEBUG(std::cout << "\t";)
  SL_DEBUG(c2->Print();)
  SL_DEBUG(std::cout << "\n";)
  SL_DEBUG(std::cout << "relative to " << m_currentEvent->get_point() << "\n";)

  CGAL_assertion(c1 != c2);
  
  Object res =
    m_traits.nearest_intersection_to_right(c1->get_last_curve(),
                                            c2->get_last_curve(), 
                                            m_currentEvent->get_point());
  if (!res.is_empty())
  {
    Point_2 xp; 
    if (!CGAL::assign(xp, res))
    {
      X_monotone_curve_2 cv;
      if (CGAL::assign(cv, res))
      {
        xp = m_traits.curve_source(cv);
        Point_2 xp1 = m_traits.curve_target(cv);
        if ( m_traits.compare_xy(xp1, xp) == LARGER )
          xp = xp1;
        SL_DEBUG(std::cout << "overlap detected\n";)
      }
    }

    SL_DEBUG(
      std::cout << " a new event is created between:\n\t";
      c1->Print();
      std::cout << "\t";
      c2->Print();
      std::cout << ")\n\t at (" 
                << xp << ")" << std::endl;
      )

    // insert the event and check if an event at this point already exists
   
   const std::pair<EventQueueIter,bool>& insert_res = 
     (m_queue->insert(EventQueueValueType(xp,0)));

    Event *e ;
    if(insert_res.second)    
    {                                   // a new event is creatd , which inidicates 
                                       // that the intersection point cannot be one 
                                       //of the end-points of two curves
      e = allocate_event(xp);
      
#ifndef NDEBUG
      e->id = m_eventId++;
#endif
      
      e->push_back_curve_to_left(c1);
      e->push_back_curve_to_left(c2);
      
      add_curve_to_right(e, c1);
      add_curve_to_right(e, c2);

      PRINT_NEW_EVENT(xp, e);
      (insert_res.first)->second = e;
 
    } 
    else   // the event already exists
    {
      SL_DEBUG(std::cout << "event already exists,updating.. (" << xp <<")\n";)
      e = (insert_res.first)->second;
    
      e->add_curve_to_left(c1);

      if ( !c1->is_end_point(e))
      { 
        add_curve_to_right(e, c1);
      }

      if(e->get_num_left_curves() == 1)
        e->push_back_curve_to_left(c2);
      else
        e->add_curve_to_left(c2); 
      
      if ( !c2->is_end_point(e) ) 
      {
        add_curve_to_right(e, c2);
      }
      SL_DEBUG(e->Print();)
    }
  } 
  // no intersetion to the right of the event
  SL_DEBUG(std::cout << "not found 2\n";)
}


//DEBUG UTILITIES
#ifndef NDEBUG
  #include<CGAL/Sweep_line_2_old/Sweep_line_2_impl_debug.C>
#endif

CGAL_END_NAMESPACE

#endif
