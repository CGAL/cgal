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
#include <algorithm>
#include <CGAL/assertions.h>
#include <CGAL/memory.h>
#include <CGAL/Sweep_line_2/Sweep_line_functors.h>
#include <CGAL/Sweep_line_2/Sweep_line_traits.h>
#include <CGAL/Sweep_line_2/Sweep_line_event.h>
#include <CGAL/Sweep_line_2/Sweep_line_subcurve.h>


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
                          std::cout << "    currentPos = " << m_sweepLinePos \
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
  typedef Point_less_functor<Point_2, Traits> PointLess;
  typedef std::map< Point_2 , Event*, PointLess> EventQueue; 
  typedef typename EventQueue::iterator EventQueueIter;
  typedef typename EventQueue::value_type EventQueueValueType;

  typedef typename Event::SubCurveIter EventCurveIter;

  typedef CurveWrap Subcurve;
  typedef std::list<Subcurve*> SubCurveList;
  typedef typename SubCurveList::iterator SubCurveListIter;
   
  typedef Status_line_curve_less_functor<Traits, Subcurve> StatusLineCurveLess;
  typedef typename StatusLineCurveLess::Compare_param CompareParams;
 
  typedef std::set<Subcurve*, StatusLineCurveLess, CGAL_ALLOCATOR(int)> StatusLine;
  typedef typename StatusLine::iterator StatusLineIter;
 
  typedef typename Allocator::template rebind<Event>      EventAlloc_rebind;
  typedef typename EventAlloc_rebind::other               EventAlloc;

  typedef typename Allocator::template rebind<Subcurve>   SubcurveAlloc_rebind;
  typedef typename SubcurveAlloc_rebind::other            SubCurveAlloc;

 


  Sweep_line_2_impl(SweepVisitor* visitor) :
      m_traits(new Traits()),
      m_sweep_line_traits(m_traits),
      m_traitsOwner(true),
      m_comp_param(new CompareParams(m_traits)),
      m_queue(new EventQueue(PointLess(m_traits))),
      m_statusLine(new StatusLine(StatusLineCurveLess(m_comp_param, &m_sweepLinePos))),
      m_xcurves(0),
      m_status_line_insert_hint(m_statusLine->begin()),
      m_num_of_subCurves(0),
#ifndef NDEBUG
      m_eventId(0),
#endif
      m_visitor(visitor)
  {
    m_visitor->attach(this);
  }


  Sweep_line_2_impl(Traits *t, SweepVisitor* visitor) :
      m_traits(t),
      m_sweep_line_traits(m_traits),
      m_traitsOwner(false),
      m_comp_param(new CompareParams(m_traits)),
      m_queue(new EventQueue(PointLess(m_traits))),
      m_statusLine(new StatusLine(StatusLineCurveLess(m_comp_param, &m_sweepLinePos))),
      m_xcurves(0),
      m_status_line_insert_hint(m_statusLine->begin()),
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
    for(CurveInputIterator iter = curves_begin; iter != curves_end; ++iter)
    {
      m_traits->curve_make_x_monotone(*iter, std::back_inserter(m_xcurves));
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
    for(CurveInputIterator iter = curves_begin; iter != curves_end; ++iter)
    {
      m_traits->curve_make_x_monotone(*iter, std::back_inserter(m_xcurves));
    }
    std::copy(xcurves_begin, xcurves_end, std::back_inserter(m_xcurves));
    m_num_of_subCurves = m_xcurves.size();
    m_subCurves = m_subCurveAlloc.allocate(m_num_of_subCurves);
    for(unsigned int index = 0; index < m_num_of_subCurves ; ++index)
    {
      init_curve(m_xcurves[index],index);
    }
  }


  
  void init_curve(X_monotone_curve_2 &curve,unsigned int index);

  /*! The main loop to calculate intersections among the curves
   *  Looping over the events in the queue, for each event we first
   *  handle the curves that are tothe left of the event point (i.e., 
   *  curves that we are done with), and then we look at the curves 
   *  to the right of the point, which means we attept to find intersections
   *  between them and their neighbours on the sweep line.
   */
  void sweep()
  {
    EventQueueIter eventIter = m_queue->begin();
    m_sweepLinePos = eventIter->first;

    while ( eventIter != m_queue->end() )
    {
      m_sweepLinePos = (eventIter->first);
      m_currentEvent = eventIter->second;

      SL_DEBUG(std::cout << "------------- " << m_sweepLinePos << " --------------"
               << std::endl;
               PrintStatusLine();
               m_currentEvent->Print(););
      
      handle_left_curves();
      m_queue->erase(eventIter);

      handle_right_curves();
      if(m_visitor->after_handle_event(m_currentEvent))
        deallocate_event(m_currentEvent);

      eventIter = m_queue->begin();
    }
  }

  protected:
  
  /*! For each left-curve, if it is the "last" subcurve, i.e., the 
   * event point is the right-edge of the original curve, the 
   * last sub curve is created and added to the result. Otherwise
   * the curve is added as is to the result.
   */
  void handle_left_curves()
  {
    m_visitor->before_handle_event(m_currentEvent);
    EventCurveIter leftCurveIter = m_currentEvent->left_curves_begin();
    const Point_2 &eventPoint = m_currentEvent->get_point();
    
    // indicates if the curve will be removed for good
    bool remove_for_good = false; 

    SL_DEBUG(std::cout << "Handling left curve" << std::endl;);
    SL_DEBUG(if( leftCurveIter != m_currentEvent->left_curves_end() )
             {
               m_currentEvent->Print();
             });

    while ( leftCurveIter != m_currentEvent->left_curves_end() )  
    {
      Subcurve *leftCurve; 
      if( (*leftCurveIter)->get_overlap_subcurve() != NULL &&
           m_traits->compare_xy(
            (*leftCurveIter)->get_overlap_subcurve()->get_left_end(),
             m_currentEvent->get_point()) == SMALLER)
      {
        leftCurve = (Subcurve*)((*leftCurveIter)->getSubcurve());
        m_currentEvent->replace_right_curve((*leftCurveIter), leftCurve);
      }
      else
      {
        leftCurve = (*leftCurveIter); 
      }

      //const X_monotone_curve_2 &cv = leftCurve->get_curve();
      //const Point_2 &lastPoint = leftCurve->get_last_point();

      if ( leftCurve->is_source(eventPoint))
      {  
        remove_for_good = true;
        /*if ( !leftCurve->is_target(lastPoint) )
        {
          X_monotone_curve_2 a,b;
          m_traits->curve_split(cv, a, b, lastPoint);
          m_visitor->add_subcurve(a, leftCurve);
        }
        else 
        {
          m_visitor->add_subcurve(cv,leftCurve);
        }*/
        m_visitor->add_subcurve(leftCurve->get_last_curve(), leftCurve);

        if(leftCurve->get_orig_subcurve1() != NULL)
        {
          leftCurve->get_orig_subcurve1()->set_overlap_subcurve(NULL);
          leftCurve->get_orig_subcurve2()->set_overlap_subcurve(NULL);

          // clip the two subcurves according to end_overlap point
          Subcurve* res;
          if(( res = (Subcurve*) (leftCurve->get_orig_subcurve1() -> clip(eventPoint))) != NULL)
          {
            add_curve_to_right(m_currentEvent, res);
            m_currentEvent->mark_internal_intersection_point();

          }

          if(( res = (Subcurve*)(leftCurve->get_orig_subcurve2() -> clip(eventPoint))) !=NULL)
          {
            add_curve_to_right(m_currentEvent, res);
            m_currentEvent->mark_internal_intersection_point();
          }
        }
      }
      else if ( leftCurve->is_target(eventPoint))
      {
        remove_for_good = true;
        /*if ( !leftCurve->is_source(lastPoint))
        {
          X_monotone_curve_2 a,b;
          m_traits->curve_split(cv, a, b, lastPoint);
          m_visitor->add_subcurve(b, leftCurve);
        } else
        {
          m_visitor->add_subcurve(cv, leftCurve);
        }*/
        m_visitor->add_subcurve(leftCurve->get_last_curve(), leftCurve);
        if(leftCurve->get_orig_subcurve1() != NULL)
        {
          leftCurve->get_orig_subcurve1()->set_overlap_subcurve(NULL);
          leftCurve->get_orig_subcurve2()->set_overlap_subcurve(NULL);

          // clip the two subcurves according to end_overlap point
          Subcurve* res;
          if(( res = (Subcurve*) (leftCurve->get_orig_subcurve1() -> clip(eventPoint)) )!= NULL)
          {
            add_curve_to_right(m_currentEvent, res);
            m_currentEvent->mark_internal_intersection_point();
          }

          if((  res =(Subcurve*)(leftCurve->get_orig_subcurve2() -> clip(eventPoint))) != NULL)
          {
            add_curve_to_right(m_currentEvent, res);
            m_currentEvent->mark_internal_intersection_point();
          }
        }

      } else { 
       /* X_monotone_curve_2 a,b;
        if ( leftCurve->is_source(lastPoint)) {
          m_traits->curve_split(leftCurve->get_last_curve(), a, b, eventPoint);
          m_visitor->add_subcurve(a, leftCurve);
        } else if ( leftCurve->is_target(lastPoint)) {
          m_traits->curve_split(leftCurve->get_last_curve(), b, a, eventPoint);
          m_visitor->add_subcurve(b, leftCurve);
        } else {
          const X_monotone_curve_2 &lastCurve = leftCurve->get_last_curve();
          if ( leftCurve->is_source_left_to_target() ) {
            m_traits->curve_split(lastCurve, a, b, eventPoint);
            m_visitor->add_subcurve(a, leftCurve);
          } else {
            m_traits->curve_split(lastCurve, b, a, eventPoint);
            m_visitor->add_subcurve(a, leftCurve);
          }*/

        X_monotone_curve_2 a,b;
        const X_monotone_curve_2 &lastCurve = leftCurve->get_last_curve();
        if ( leftCurve->is_source_left_to_target() ) 
        {
          m_traits->curve_split(lastCurve, a, b, eventPoint);
          m_visitor->add_subcurve(a, leftCurve);
        }
        else
        {
          m_traits->curve_split(lastCurve, b, a, eventPoint);
          m_visitor->add_subcurve(a, leftCurve);
        }
        leftCurve->set_last_point(eventPoint);
        leftCurve->set_last_curve(b); 
      }
      
      // remove curve from the status line (also checks intersection 
      // between the neighbouring curves,only if the curve is removed for good)
      remove_curve_from_status_line(leftCurve, remove_for_good);

      ++leftCurveIter;
    }
    SL_DEBUG(std::cout << "Handling left curve END" << std::endl;)
      
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
    int numRightCurves = m_currentEvent->get_num_right_curves();
    if ( numRightCurves == 0 )
      return;
      
    
    int numLeftCurves = m_currentEvent->get_num_left_curves();
    /* this block takes care of 
    //
    //           /
    //          /
    //       --------
    //          \
    //           \
    */
    if ( numLeftCurves == 0 ) 
    {
      SL_DEBUG(std::cout << " - handling special case " << std::endl;)
                                                     
      EventCurveIter currentOne = m_currentEvent->right_curves_begin();
      StatusLineIter slIter = m_statusLine->lower_bound(*currentOne);
      if ( slIter != m_statusLine->end() &&
           CurveStartsAtCurve(*currentOne, *slIter))
      {
        Subcurve* sc = *slIter;
        X_monotone_curve_2 last_curve = sc->get_last_curve();
        m_currentEvent->add_curve_to_left(sc); 
        bool is_overlap = (add_curve_to_right(m_currentEvent, sc)).first;
        m_currentEvent->mark_internal_intersection_point(); 
        X_monotone_curve_2 a,b;
        if ( sc->is_source_left_to_target() )
        {
          m_traits->curve_split(last_curve, a, b,m_currentEvent->get_point());
        }
        else
        {
          m_traits->curve_split(last_curve, b, a, 
            m_currentEvent->get_point());
        }
        if(is_overlap)
        {
          m_status_line_insert_hint = slIter; ++m_status_line_insert_hint; 
          m_statusLine->erase(slIter);
        }
        else
        {
          sc->set_last_point(m_currentEvent->get_point());
          sc->set_last_curve(b);
        }
        m_visitor->add_subcurve(a, sc);
      }
      else
      {
        if( slIter != m_statusLine->begin() &&
          CurveStartsAtCurve(*currentOne, *(--slIter)))
        {
          Subcurve* sc = *slIter;
          X_monotone_curve_2 last_curve = sc->get_last_curve();
          m_currentEvent->add_curve_to_left(sc); 
          bool is_overlap =  (add_curve_to_right(m_currentEvent,sc)).first;
          m_currentEvent->mark_internal_intersection_point(); 
          X_monotone_curve_2 a,b;
          if ( sc->is_source_left_to_target() ) 
          {
            m_traits->curve_split(last_curve, a, b, 
              m_currentEvent->get_point());
          }
          else
          {
            m_traits->curve_split(last_curve, b, a, 
              m_currentEvent->get_point());
          }
          if(is_overlap)
          {
            m_status_line_insert_hint = slIter; ++m_status_line_insert_hint; 
            m_statusLine->erase(slIter);
          }
          else
          {
            ++numRightCurves;
            sc->set_last_point(m_currentEvent->get_point());
            sc->set_last_curve(b);
          }
          m_visitor->add_subcurve(a,sc);
        }
      }
    }
    if(numRightCurves == 1)
    {
      SL_DEBUG(std::cout << " - beginning of curve " << std::endl;);
      SL_DEBUG(
        Subcurve *tmp1 = *(m_currentEvent->right_curves_begin());
        PRINT_INSERT(tmp1);
               );

      StatusLineIter slIter = m_statusLine->insert(m_status_line_insert_hint, 
                                *(m_currentEvent->right_curves_begin()));
      
      (*(m_currentEvent->right_curves_begin()))->set_hint(slIter); //xxx
      m_status_line_insert_hint = slIter; ++m_status_line_insert_hint;
     
      SL_DEBUG(PrintStatusLine(););
     
      // if this is the only curve on the status line, nothing else to do
      if ( m_statusLine->size() == 1 )
        return;

      StatusLineIter prev = slIter;
      StatusLineIter next = slIter;
      ++next;
      if ( slIter != m_statusLine->begin() )
      {
        --prev;
        intersect(*prev, *(m_currentEvent->right_curves_begin()));
      }
      if ( next != m_statusLine->end() )
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

      StatusLineIter slIter = m_statusLine->insert(m_status_line_insert_hint, 
                                                    *firstOne);
      ((Subcurve*)(*firstOne))->set_hint(slIter);
        
      SL_DEBUG(PrintStatusLine(););
      if ( slIter != m_statusLine->begin() )
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
        ++slIter;
        slIter = m_statusLine->insert(slIter, *currentOne);
        ((Subcurve*)(*currentOne))->set_hint(slIter);
          
        SL_DEBUG(PrintStatusLine(););
      
        intersect(*prevOne, *currentOne);
        prevOne = currentOne;
        ++currentOne;
      }
      m_status_line_insert_hint = slIter; ++m_status_line_insert_hint;
        
      lastOne = currentOne; --lastOne;
        
      SL_DEBUG(PrintStatusLine(););
      StatusLineIter next = slIter; ++next;
      if ( next != m_statusLine->end() )
        intersect( *prevOne,*next);
    }
  }



  std::pair<bool,SubCurveListIter> add_curve_to_right(Event* event, Subcurve* curve)
  {
    std::pair<bool, SubCurveListIter> pair_res = event->add_curve_to_right(curve);
    if(pair_res.first == true) //overlap
    {
      SubCurveListIter iter = pair_res.second;
      
       X_monotone_curve_2 overlap_cv;
       Object cv_obj =
        m_traits->nearest_intersection_to_right(curve->get_curve(),
                                                (*iter)->get_curve(), 
                                                event->get_point());
      if (CGAL::assign(overlap_cv, cv_obj))
      {
         //alocate new Subcure for the overlap
         Subcurve *overlap_sc = m_subCurveAlloc.allocate(1);
         m_subCurveAlloc.construct(overlap_sc,m_masterSubcurve);
         overlap_sc->init(overlap_cv);
         m_overlap_subCurves.push_back(overlap_sc);
         
         //get the right end of overlap_cv
         Point_2 end_overlap = m_traits->curve_target(overlap_cv);
         if(m_traits->compare_xy(end_overlap,
                                 m_traits->curve_source(overlap_cv)) == SMALLER)
         {
           end_overlap = m_traits->curve_source(overlap_cv);
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
    return pair_res;
  }

  


  // utility methods 
  bool intersect(Subcurve *c1, Subcurve *c2);
  void remove_curve_from_status_line(Subcurve *leftCurve,bool remove_for_good);

 
  
#ifndef NDEBUG
  void PrintEventQueue();
  void PrintSubCurves();
  void PrintStatusLine();
#endif

protected:

  /*! a pointer to a traits object */
  Traits *m_traits;

  /*! an object that holds a static trait object 
   *  to be used by events and subcurves
   */
  Sweep_line_traits<Traits> m_sweep_line_traits;

  /*! an indication to whether the traits should be deleted in the destructor
   */
  bool m_traitsOwner;

  /*! The current position (in X) of the sweep line */
  Point_2 m_sweepLinePos;

  /*! a struct that holds params associated with the curve compare functor */
  CompareParams *m_comp_param;

  /*! the queue of events (intersection points) to handle */
  EventQueue *m_queue;

  /*! The subcurves array */
  Subcurve *m_subCurves;

  /*! The status line (Y-str) */
  StatusLine *m_statusLine;
 
  /*! if non x-monotone are specified, this hold the x-monotone 
    curves created when splitting them into x-monotone curves. */
  std::vector<X_monotone_curve_2> m_xcurves;

  /*! a pointer to the current event */
  Event *m_currentEvent;

  /*! An iterator of the  status line that is used as a hint for inserts. */
  StatusLineIter m_status_line_insert_hint;

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

 

  bool CurveStartsAtCurve(Subcurve *one, Subcurve *two)
  {
    SL_DEBUG(std::cout << "CurveStartsAtCurve: \n";);
    SL_DEBUG(std::cout << one->get_curve() << "\n" << two->get_curve() 
                       << "\n";);

    if ( m_traits->curve_compare_y_at_x(one->get_left_end(), 
                                        two->get_curve()) == EQUAL )
      return true;
    return false;
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
  if ( m_traitsOwner ) delete m_traits;
  
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
  delete m_statusLine;
  delete m_comp_param;
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
 
  const Point_2 &left_end  =  (m_subCurves+j)->get_left_end();
  const Point_2 &right_end =  (m_subCurves+j)->get_right_end();
  

    //handle the right end of the curve
    const std::pair<EventQueueIter, bool>& insertion_res =
      m_queue->insert(EventQueueValueType(right_end,0));  //the insertion return value

    if(insertion_res.second == true)  // its a new event
    {
      e =  m_eventAlloc.allocate(1); 
      m_eventAlloc.construct(e, m_masterEvent);
      e->init(right_end );
     

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
      e =  m_eventAlloc.allocate(1);
      m_eventAlloc.construct(e, m_masterEvent);
      e -> init(left_end);
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
    m_statusLine->erase(sliter);
    SL_DEBUG(std::cout << "remove_curve_from_status_line Done\n";)
    return;
  }

  CGAL_assertion(sliter!=m_statusLine->end());
  StatusLineIter end = m_statusLine->end(); --end;
  if ( sliter != m_statusLine->begin() && sliter != end ) 
  {
    StatusLineIter prev = sliter; --prev;
    StatusLineIter next = sliter; ++next;
    
    // intersect *next with  *prev 
    intersect(*prev, *next);
  }
  m_statusLine->erase(sliter);
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
inline bool 
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


  if ( c1 == c2 )
  {
    SL_DEBUG(std::cout << "same curve, returning....\n";)
    return false;
  }
 
  const X_monotone_curve_2 &cv1 = c1->get_curve();
  const X_monotone_curve_2 &cv2 = c2->get_curve();

  bool isOverlap = false;

  Object res =
    m_traits->nearest_intersection_to_right(cv1, cv2, 
                                            m_currentEvent->get_point());
  if (!res.is_empty())
  {
    Point_2 xp; 
    if (!CGAL::assign(xp, res))
    {
      X_monotone_curve_2 cv;
      if (CGAL::assign(cv, res))
      {
        xp = m_traits->curve_source(cv);
        Point_2 xp1 = m_traits->curve_target(cv);
        if ( m_traits->compare_xy(xp1, xp) == LARGER )
          xp = xp1;
        SL_DEBUG(std::cout << "overlap detected\n";)
        isOverlap = true;
      }
    }

    SL_DEBUG(
      std::cout << " a new event is created between:\n\t";
      c1->Print();
      std::cout << "\t";
      c2->Print();
      std::cout << "\trelative to ("
                << m_sweepLinePos << ")\n\t at (" 
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
      e =  m_eventAlloc.allocate(1);
      m_eventAlloc.construct(e, m_masterEvent);
      e -> init(xp);
      
#ifndef NDEBUG
      e->id = m_eventId++;
#endif
      
      e->add_curve_to_left(c1);
      if(e->get_num_left_curves() == 1)
        e->push_back_curve_to_left(c2);
      else
        e->add_curve_to_left(c2); 

      add_curve_to_right(e, c1);
      add_curve_to_right(e, c2);

      e->mark_internal_intersection_point();       
      PRINT_NEW_EVENT(xp, e);
      (insert_res.first)->second = e;
 
      return isOverlap;
    } 
    else   // the event already exists
    {
      SL_DEBUG(std::cout << "event already exists,updating.. (" << xp <<")\n";)
      e = (insert_res.first)->second;
    
      e->add_curve_to_left(c1);
      if ( !c1->is_end_point(xp)) 
      { 
        add_curve_to_right(e, c1);
        e->mark_internal_intersection_point();
      }

      if(e->get_num_left_curves() == 1)
        e->push_back_curve_to_left(c2);
      else
        e->add_curve_to_left(c2); 
      
      if ( !c2->is_end_point(xp) ) 
      {
        add_curve_to_right(e, c2);
        e->mark_internal_intersection_point();
      }
      SL_DEBUG(e->Print();)
    }
    return isOverlap;
  } 
  // no intersetion to the right of the event
  SL_DEBUG(std::cout << "not found 2\n";)
  return isOverlap;
}




//DEBUG UTILITIES
#ifndef NDEBUG
  #include<CGAL/Sweep_line_2/Sweep_line_2_impl_debug.C>
#endif



CGAL_END_NAMESPACE

#endif
