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
#ifndef CGAL_SWEEP_LINE_2_IMPL_H
#define CGAL_SWEEP_LINE_2_IMPL_H


#include <CGAL/Sweep_line_2/Sweep_line_functors.h>
#include <CGAL/assertions.h>
#include <map>
#include <set>
#include <CGAL/memory.h>
#include <CGAL/Sweep_line_2/Sweep_line_traits.h>

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

#include <list>

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

  There are two main functionalities supported by this algorithm:
  1. calculate the non intersecting curves that are a product of 
     intersecting a set of input curves
  2. calculate all the intersection points between the curves specified

  An extension of this algorithm is used to produce a planar map containing
  the input curves efficiently. \sa Pmwx_aggregate_insert.

  General flow:
  After the initialization stage, the events are handled from left to right.

  For each event

    First pass - handles special cases in which curves start or end 
                 at the interior of another curve
    Handle left curves - iterate over the curves that intersect 
                 at the event point and defined to the left of the 
                 event. This is mostly where the output (sub curves 
                 or intersection points) is produced.
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
    out - always refer to the output iterator

*/

template <class CurveInputIterator,  class SweepLineTraits_2,
          class SweepEvent, class CurveWrap,
          typename EventAlloc    = CGAL_ALLOCATOR (SweepEvent),
          typename SubCurveAlloc = CGAL_ALLOCATOR (CurveWrap) >
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
 
  class  SweepLineGetSubCurves {};
  class  SweepLineGetPoints {};
  class  SweepLineGetInterCurveList {};
  class  SweepLinePlanarmap {};

  Sweep_line_2_impl()  : m_traits(new Traits()),
                         m_sweep_line_traits(m_traits),
                         m_traitsOwner(true),
                         m_includeEndPoints(true),
                         m_queue(NULL),
                         m_statusLine(NULL),
                         m_comp_param(NULL),
                         m_xcurves(0),
                         m_found_intersection(false),
                         is_first_point(true),
                         m_num_of_subCurves(0)
  {
  }


  Sweep_line_2_impl(Traits *t) : m_traits(t),
                                 m_sweep_line_traits(m_traits),
                                 m_traitsOwner(false),
                                 m_includeEndPoints(true),
                                 m_queue(NULL),
                                 m_statusLine(NULL),
                                 m_comp_param(NULL),
                                 m_xcurves(0),
                                 m_found_intersection(false),
                                 is_first_point(true),
                                 m_num_of_subCurves(0)
  {
  }

  virtual ~Sweep_line_2_impl();

  /*!
   *  Given a container of curves, this function returns a list of curves
   *  that are created by intersecting the input curves.
   *  \param curves_begin the input iterator that points to the first curve 
   *                      in the range.
   *  \param curves_end the input past-the-end iterator of the range.
   *  \param subcurves an iterator to the first curve in the range
   *                   of curves created by intersecting the input curves.
   *  \param overlapping indicates whether overlapping curves should be 
   *                   reported once or multiple times. If false, the 
   *                   overlapping curves are reported once only.
   */
  template <class OutpoutIterator>
  void  get_subcurves(CurveInputIterator begin, CurveInputIterator end, 
                      OutpoutIterator subcurves, bool overlapping = false)
  { 
    init(begin, end);
    SL_DEBUG(
      PrintSubCurves();
      PrintEventQueue();
    )
    m_overlapping = overlapping;
    sweep(subcurves, SweepLineGetSubCurves());
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
    init(begin, end);
    SL_DEBUG(
      PrintSubCurves();
      PrintEventQueue();
    )
    m_includeEndPoints = includeEndPoints;
    m_found_intersection = false;
    sweep(points, SweepLineGetPoints());
  }

 /*!
  *  Given a range of curves, this function returns an iterator 
  *  to the beginning of a range that contains the list of curves 
  *  for each intersection point between any two curves in the 
  *  specified range.
  *  The intersections are calculated using the sweep algorithm.
  *  \param begin the input iterator that points to the first curve 
  *                      in the range.
  *  \param end the input past-the-end iterator of the range.
  *  \param intersecting_curves an iterator to the output
  *  \param endpoints if true, the end points of the curves are reported
  *                   as intersection points. Defaults to true.
  */
  template <class OutputIterator>
  void  get_intersecting_curves(CurveInputIterator begin, 
                                CurveInputIterator end, 
                                OutputIterator intersecting_curves,
                                bool endpoints = true)
  { 
    // TODO - implement...
    SweepLineGetInterCurveList tag;
    (void) tag;
  }

  bool do_curves_intersect(CurveInputIterator begin, 
                            CurveInputIterator end)
  {
    init(begin, end);
    SL_DEBUG(
      PrintSubCurves();
      PrintEventQueue();
    )
    m_includeEndPoints = false;
    std::vector<Point_2> dummy;
    m_stop_at_first_int = true;
    m_found_intersection = false;
    sweep(std::back_inserter(dummy), SweepLineGetPoints(), true);
    return m_found_intersection;
  }


  static Traits* get_traits()
  {
    static Traits* s_traits = 0;
    if(!s_traits)
      return new Traits;
    return s_traits;
  }

protected:


  void init(CurveInputIterator begin, CurveInputIterator end);
  void init_curve(X_monotone_curve_2 &curve,unsigned int index);

  /*! The main loop to calculate intersections among the curves
   *  Looping over the events in the queue, for each event we first
   *  handle the curves that are tothe left of the event point (i.e., 
   *  curves that we are done with), and then we look at the curves 
   *  to the right of the point, which means we attept to find intersections
   *  between them and their neighbours on the sweep line.
   */

  template <class OutpoutIterator, class Op>
  void sweep(OutpoutIterator out, Op tag, bool stop_at_first_int=false)
  {
    EventQueueIter eventIter = m_queue->begin();
    m_sweepLinePos = eventIter->first;

    while ( eventIter != m_queue->end() )
    {
      const Point_2 &p = (eventIter->first);
      m_sweepLinePos = p;
   
      m_currentEvent = eventIter->second;
      SL_DEBUG(std::cout << "------------- " << p << " --------------"
               << std::endl;
               PrintStatusLine();
               m_currentEvent->Print();
      )
      
     
      handle_left_curves(out, tag); 
      m_queue->erase(eventIter);
      handle_right_curves(out, tag);

      m_eventAlloc.destroy(m_currentEvent);
      m_eventAlloc.deallocate(m_currentEvent,1);

      eventIter = m_queue->begin();
    }

    if ( stop_at_first_int && m_found_intersection )
      return;
  }

  
  /*! For each left-curve, if it is the "last" subcurve, i.e., the 
   * event point is the right-edge of the original curve, the 
   * last sub curve is created and added to the result. Otherwise
   * the curve is added as is to the result.
   */
  template <class OutpoutIterator>
  void handle_left_curves(OutpoutIterator out, 
                        SweepLineGetSubCurves &tag)
  {
    SL_DEBUG(std::cout << "Handling left curve" << std::endl;)
    SL_DEBUG(m_currentEvent->Print();)

    EventCurveIter leftCurveIter = m_currentEvent->left_curves_begin();
    const Point_2 &eventPoint = m_currentEvent->get_point();

    while ( leftCurveIter != m_currentEvent->left_curves_end() )  
    {
      Subcurve *leftCurve = *leftCurveIter; 
      const X_monotone_curve_2 &cv = leftCurve->get_curve();
      const Point_2 &lastPoint = leftCurve->get_last_point();

      if ( leftCurve->is_source(eventPoint))
      {      
        if ( !leftCurve->is_target(lastPoint) )
        {
          X_monotone_curve_2 a,b;
          m_traits->curve_split(cv, a, b, lastPoint);
          add_curve_to_output(a, leftCurve, out);
        } else {
          add_curve_to_output(cv, leftCurve, out);
        }
      } else if ( leftCurve->is_target(eventPoint))
      {
        if ( !leftCurve->is_source(lastPoint))
        {
          X_monotone_curve_2 a,b;
          m_traits->curve_split(cv, a, b, lastPoint);
          add_curve_to_output(b, leftCurve, out);
        } else {
          add_curve_to_output(cv, leftCurve, out);
        }

      } else { 
        X_monotone_curve_2 a,b;
        if ( leftCurve->is_source(lastPoint)) {
          m_traits->curve_split(cv, a, b, eventPoint);
          add_curve_to_output(a, leftCurve, out);
        } else if ( leftCurve->is_target(lastPoint)) {
          m_traits->curve_split(cv, b, a, eventPoint);
          add_curve_to_output(a, leftCurve, out);
        } else {
          const X_monotone_curve_2 &lastCurve = leftCurve->get_last_curve();
          if ( leftCurve->is_source_left_to_target() ) {
            m_traits->curve_split(lastCurve, a, b, eventPoint);
            add_curve_to_output(a, leftCurve, out);
          } else {
            m_traits->curve_split(lastCurve, b, a, eventPoint);
            add_curve_to_output(a, leftCurve, out);
          }
        }
        leftCurve->set_last_point(eventPoint);
        leftCurve->set_last_curve(b); 
      
      }
      
      // remove curve from the status line (also checks intersection 
      // between the neighbouring curves)
      remove_curve_from_status_line(leftCurve);

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
  template <class OutpoutIterator>
  void handle_right_curves(OutpoutIterator out, SweepLineGetSubCurves &tag)
  {
    SL_DEBUG(std::cout << "Handling right curves (" ;)
    SL_DEBUG(std::cout << m_currentEvent->get_point() << ")\n";)
    int numRightCurves = m_currentEvent->get_num_right_curves();
    if ( numRightCurves == 0 )
      return;
      
    if ( numRightCurves == 1 )
    {
      SL_DEBUG(std::cout << " - beginning of curve " << std::endl;)
          
      SL_DEBUG(
            Subcurve *tmp1 = *(m_currentEvent->right_curves_begin());
            PRINT_INSERT(tmp1);
            )
          
      StatusLineIter slIter = m_statusLine->insert(
                                m_status_line_insert_hint, 
                                *(m_currentEvent->right_curves_begin()));
      
      (*(m_currentEvent->right_curves_begin()))->set_hint(slIter); //xxx
      m_status_line_insert_hint = slIter; ++m_status_line_insert_hint;

      SL_DEBUG(PrintStatusLine();)

      // if this is the only curve on the status line, nothing else to do
      if ( m_statusLine->size() == 1 )
        return;

      StatusLineIter prev = slIter;
      StatusLineIter next = slIter;
      ++next;

      SubCurveList mylist;
      if ( slIter != m_statusLine->begin() )
      {
        --prev;
        
        StatusLineIter tmp = prev;
        mylist.push_back(*prev);
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
        mylist.push_back(*next);
        ++tmp;
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
      intersect_curve_group(*(m_currentEvent->right_curves_begin()), 
                            mylist, out);
      
    } else
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
      if ( numLeftCurves == 0 ) {
        
        SL_DEBUG(std::cout << " - handling special case " << std::endl;)
        StatusLineIter slIter;
                                                     
        EventCurveIter currentOne = m_currentEvent->right_curves_begin();
        while ( currentOne != m_currentEvent->right_curves_end() ) {
          slIter = m_statusLine->lower_bound(*currentOne);
          if ( slIter != m_statusLine->end() ) {
            Subcurve *c = *slIter;
            if ( CurveStartsAtCurve(*currentOne, c)) {
              m_currentEvent->add_curve_to_left(c); 
              m_currentEvent->add_curve_to_right(c);
        m_currentEvent->mark_internal_intersection_point(); 
              X_monotone_curve_2 a,b;
              if ( c->is_source_left_to_target() ) {
                m_traits->curve_split(c->get_last_curve(), a, b, 
                                      m_currentEvent->get_point());
              } else {
                m_traits->curve_split(c->get_last_curve(), b, a, 
                                      m_currentEvent->get_point());
              }
              c->set_last_point(m_currentEvent->get_point());
              c->set_last_curve(b);
              
              add_curve_to_output(a, c, out);
              break;
            }
          }
          currentOne++;
        }
      }
      // end block ...
      
      
      
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
      (*firstOne)->set_hint(slIter);
      
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
        
        intersect_curve_group(*slIter, prevlist, out);
      }
      currentlist.push_back(*firstOne);
      
      EventCurveIter currentOne = firstOne; ++currentOne;
      EventCurveIter prevOne = firstOne;
      
      while ( currentOne != rightCurveEnd )
      {
        PRINT_INSERT(*currentOne);
        ++slIter;
        slIter = m_statusLine->insert(slIter, *currentOne);
        (*currentOne)->set_hint(slIter);
        
        SL_DEBUG(PrintStatusLine(););
        if ( do_curves_overlap(*currentOne, *prevOne))
        {
          intersect_curve_group(*currentOne, currentlist, out);
          currentlist.push_back(*currentOne);
        } else {
          prevlist = currentlist;
          currentlist.clear();
          currentlist.push_back(*currentOne);
        }
        
        intersect_curve_group(*currentOne, prevlist, out);
        prevOne = currentOne;
        ++currentOne;
      }
      m_status_line_insert_hint = slIter; ++m_status_line_insert_hint;
      
      lastOne = currentOne; --lastOne;
      
      SL_DEBUG(PrintStatusLine();)
        StatusLineIter next = slIter; ++next;
      if ( next != m_statusLine->end() ) {
        intersect_curve_group(*next, currentlist, out, true);
        StatusLineIter tmp = next; ++tmp;
        while ( tmp != m_statusLine->end() ) 
        {
          if ( do_curves_overlap(*next, *tmp))
          {
            intersect_curve_group(*tmp, currentlist, out, true);
            ++tmp;
          }
          else
            break;
        }
      }
    }
  }
  

  template <class OutpoutIterator>
  void handle_right_curves(OutpoutIterator out, SweepLineGetPoints &tag)
  {
    SL_DEBUG(std::cout << "Handling right curves (" ;);
    SL_DEBUG(std::cout << m_currentEvent->get_point() << ")\n";);
    int numRightCurves = m_currentEvent->get_num_right_curves();
    if ( numRightCurves == 0 )
      return;
    
    if ( numRightCurves == 1 )
    {
      SL_DEBUG(std::cout << " - beginning of curve " << std::endl;);

      SL_DEBUG(
        Subcurve *tmp1 = *(m_currentEvent->right_curves_begin());
        PRINT_INSERT(tmp1);
        );

      StatusLineIter slIter = m_statusLine->insert(
                                      m_status_line_insert_hint, 
                                      *(m_currentEvent->right_curves_begin()));
      
      (*(m_currentEvent->right_curves_begin()))->set_hint(slIter); 
      m_status_line_insert_hint = slIter; ++m_status_line_insert_hint;
      
      SL_DEBUG(PrintStatusLine(););
      
      // if this is the only curve on the status line, nothing else to do
      if ( m_statusLine->size() == 1 )
        return;
      
      StatusLineIter prev = slIter;
      StatusLineIter next = slIter;
      ++next;

      SubCurveList mylist;
      if ( slIter != m_statusLine->begin() )
      {
        --prev;
        
        if ( CurveStartsAtCurve(*slIter, *prev) && !m_includeEndPoints) {
          SL_DEBUG(std::cout << "Reporting point (4): " 
                             << (*slIter)->get_left_end() 
                             << "\n";);
          add_point_to_output((*slIter)->get_left_end(), out);
          m_found_intersection = true;
        }
        StatusLineIter tmp = prev;
        mylist.push_back(*prev);
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
        if ( CurveStartsAtCurve(*slIter, *next)  && !m_includeEndPoints) {
          SL_DEBUG(std::cout << "Reporting point (5): " 
                   << (*slIter)->get_left_end() 
                   << "\n";);
          add_point_to_output((*slIter)->get_left_end(), out);
          m_found_intersection = true;
        }
        StatusLineIter tmp = next;
        mylist.push_back(*next);
        ++tmp;
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
      intersect_curve_group(*(m_currentEvent->right_curves_begin()), mylist);
      
    } else
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
      if ( numLeftCurves == 0 ) {
        
        SL_DEBUG(std::cout << " - handling special case " << std::endl;);
        StatusLineIter slIter;
        EventCurveIter currentOne = m_currentEvent->right_curves_begin();
        while ( currentOne != m_currentEvent->right_curves_end() ) {
          slIter = m_statusLine->lower_bound(*currentOne);
          if ( slIter != m_statusLine->end() ) {
            Subcurve *c = *slIter;
            if ( CurveStartsAtCurve(*currentOne, c) && !m_includeEndPoints) {
              SL_DEBUG(std::cout << "Reporting point (8): "
                                 << (*currentOne)->get_left_end()
                                 << "\n";);
              add_point_to_output((*currentOne)->get_left_end(), out);
              m_found_intersection = true;
              break;
            }
          }
          currentOne++;
        }
      }
      // end block ..
      
      
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
      (*firstOne)->set_hint(slIter);
      
      SL_DEBUG(PrintStatusLine(););
      if ( slIter != m_statusLine->begin() )
      { 
        StatusLineIter prev = slIter; --prev;
        
        if ( CurveStartsAtCurve(*slIter, *prev) && !m_includeEndPoints) {
          SL_DEBUG(std::cout << "Reporting point (6): " 
                   << (*slIter)->get_left_end() 
                   << "\n";);
          add_point_to_output((*slIter)->get_left_end(), out);
          m_found_intersection = true;
        }
        
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
        
        intersect_curve_group(*slIter, prevlist);
      }
      currentlist.push_back(*firstOne);
      
      EventCurveIter currentOne = firstOne; ++currentOne;
      EventCurveIter prevOne = firstOne;
      
      while ( currentOne != rightCurveEnd )
      {
        PRINT_INSERT(*currentOne);
        ++slIter;
        slIter = m_statusLine->insert(slIter, *currentOne);
        (*currentOne)->set_hint(slIter);
        
        SL_DEBUG(PrintStatusLine(););
        if ( do_curves_overlap(*currentOne, *prevOne))
        {
          intersect_curve_group(*currentOne, currentlist);
          currentlist.push_back(*currentOne);
        } else {
          prevlist = currentlist;
          currentlist.clear();
          currentlist.push_back(*currentOne);
        }
        
        intersect_curve_group(*currentOne, prevlist);
        prevOne = currentOne;
        ++currentOne;
      }
      m_status_line_insert_hint = slIter; ++m_status_line_insert_hint;
      
      lastOne = currentOne; --lastOne;
      PRINT_INSERT(*lastOne);
      
      SL_DEBUG(PrintStatusLine(););
      StatusLineIter next = slIter; ++next;
      if ( next != m_statusLine->end() ) {
        
        if ( CurveStartsAtCurve(*slIter, *next)  && !m_includeEndPoints) {
          SL_DEBUG(std::cout << "Reporting point (7): " 
                             << (*slIter)->get_left_end() 
                               << "\n";);
          add_point_to_output((*slIter)->get_left_end(), out);
          m_found_intersection = true;
        }
        
        intersect_curve_group(*next, currentlist);
        StatusLineIter tmp = next; ++tmp;
        while ( tmp != m_statusLine->end() ) 
        {
          if ( do_curves_overlap(*next, *tmp))
          {
            intersect_curve_group(*tmp, currentlist);
            ++tmp;
          }
          else
            break;
        }
      }
    }
  }
  


  // utility methods 
  bool intersect(Subcurve *c1, Subcurve *c2);
  void intersect_curve_group(Subcurve *c1, SubCurveList &mylist);


  // reverse = false ==> we check if the curve starts at the list of curves
  // reverse = true ==> we check if any of the curves in the list start at 
  // the single curve
  template <class OutpoutIterator>
  void intersect_curve_group(Subcurve *c1, SubCurveList &mylist,
                           OutpoutIterator out, bool reverse=false)
  {
    m_tmpOut.clear();
    SL_DEBUG(std::cout << "intersect_curve_group (with out)\n";)
    SL_DEBUG(std::cout << "intersecting with " << mylist.size()
             << " curves\n";)
    SubCurveListIter i = mylist.begin();
    while ( i != mylist.end())
    {
      bool flag;
      if ( reverse ) {
        flag = CurveStartsAtCurve(*i, c1);
        if ( flag && (c1->get_last_point() != m_currentEvent->get_point()) ) {
          SL_DEBUG(std::cout << "CurveStartsAtCurve 3 \n";)
          m_currentEvent->add_curve_to_right(c1);
          m_currentEvent->add_curve_to_left(c1);
    m_currentEvent->mark_internal_intersection_point(); 
          X_monotone_curve_2 a,b;
          SL_DEBUG(std::cout << "splitting " << (c1)->get_last_curve() 
                             << " at " 
                             << m_currentEvent->get_point() << "\n";)
          if ( (c1)->is_source_left_to_target() ) 
            m_traits->curve_split((c1)->get_last_curve(), a, b, 
                                  m_currentEvent->get_point());
          else
            m_traits->curve_split((c1)->get_last_curve(), b, a, 
                                  m_currentEvent->get_point());
          (c1)->set_last_point(m_currentEvent->get_point());
          (c1)->set_last_curve(b); 
          (c1)->set_last_subcurve(a); 
          m_tmpOut.push_back(c1);
        }
      }
      else
      {
        flag = CurveStartsAtCurve(c1, *i);
        if ( flag && ((*i)->get_last_point() != m_currentEvent->get_point())) {
          
          SL_DEBUG(std::cout << "CurveStartsAtCurve 3 \n";)
          m_currentEvent->add_curve_to_right(*i);
          m_currentEvent->add_curve_to_left(*i);
    m_currentEvent->mark_internal_intersection_point(); 
          X_monotone_curve_2 a,b;
          SL_DEBUG(std::cout << "splitting " << (*i)->get_last_curve() 
                             << " at " 
                             << m_currentEvent->get_point() << "\n";)
          if ( (*i)->is_source_left_to_target() ) 
            m_traits->curve_split((*i)->get_last_curve(), a, b, 
                                  m_currentEvent->get_point());
          else
            m_traits->curve_split((*i)->get_last_curve(), b, a, 
                                  m_currentEvent->get_point());
          (*i)->set_last_point(m_currentEvent->get_point());
          (*i)->set_last_curve(b); 
          (*i)->set_last_subcurve(a); 
          m_tmpOut.push_back(*i);
          
        }
      }
      
      intersect(c1, *i);
      ++i;
    }
    
    for ( SubCurveListIter iter = m_tmpOut.begin() ; iter != m_tmpOut.end();
          ++iter) {
      add_curve_to_output((*iter)->get_last_subcurve(), *iter, out);
    }
  }

  void remove_curve_from_status_line(Subcurve *leftCurve);

 
  bool do_curves_overlap(Subcurve *c1, Subcurve *c2);
  bool similar_curves(const X_monotone_curve_2 &a, 
                      const X_monotone_curve_2 &b);
 

  /*! Adds a new curve to the output list. If the overlapping flag is false,
      each unique curve is reported only once.

      @param out an iterator to the output container
      @param cv the curve to be added
  */
  template <class OutpoutIterator>
  void add_curve_to_output(const X_monotone_curve_2 &cv, Subcurve *curve, 
                           OutpoutIterator out)
  {
    static Subcurve *prevCurve = 0;
    static X_monotone_curve_2 prevXCv;

    if ( m_overlapping ) {
      *out = cv;
      ++out;
    } else {
      if ( prevCurve && similar_curves(cv, prevXCv)) {
        SL_DEBUG(std::cout << " curve already reported... " << std::endl;)
        return;
      }
      prevCurve = curve;
      prevXCv = cv;
      *out = cv;
      ++out;
    }
  }
  template <class OutpoutIterator>
  void add_point_to_output(const Point_2 p, OutpoutIterator out) {

    SL_DEBUG(std::cout << "Reporting point " << p;)
    if ( is_first_point ) {
      is_first_point = false;
      *out = p;
      m_lastReportedPoint = p;
      ++out;
      SL_DEBUG(std::cout << " YES - first time\n";)
    } else {
      if ( m_lastReportedPoint != p ) {
        *out = p;
        ++out;
        m_lastReportedPoint = p;
        SL_DEBUG(std::cout << " YES\n";)
      } else {
        SL_DEBUG(std::cout << " NO\n";)
      }
    }
  }


 

  /*! For each left-curve, if it is the "last" subcurve, i.e., the 
   * event point is the right-edge of the original curve, the 
   * last sub curve is created and added to the result. Otherwise
   * the curve is added as is to the result.
   */
  template <class OutpoutIterator>
  void handle_left_curves(OutpoutIterator out, SweepLineGetPoints &tag)
  {
    SL_DEBUG(std::cout << "Handling left curve" << std::endl;)
    SL_DEBUG(m_currentEvent->Print();)
    const Point_2 &eventPoint = m_currentEvent->get_point();
    if ( !m_currentEvent->has_left_curves())
    {
      if (m_includeEndPoints || 
          m_currentEvent->is_internal_intersection_point())
      {
        SL_DEBUG(std::cout << "Reporting point (2): " << eventPoint << "\n";)
          //*out = eventPoint; ++out;    
        add_point_to_output(eventPoint, out);
        m_found_intersection = true;
      }
      return;
    }

    // delete the curve from the status line
    EventCurveIter leftCurveIter = m_currentEvent->left_curves_begin();

    while ( leftCurveIter != m_currentEvent->left_curves_end() )
    {
      // before deleting check new neighbors that will become after deletion
      remove_curve_from_status_line(*leftCurveIter);
      PRINT_ERASE((*leftCurveIter));
      Subcurve *leftCurve = *leftCurveIter; 
      leftCurve->set_last_point(eventPoint);
      ++leftCurveIter;
    }

    if ( m_includeEndPoints || 
         m_currentEvent->is_internal_intersection_point() )
    {        
      SL_DEBUG(std::cout << "Reporting point (3): " << eventPoint << "\n";)
      add_point_to_output(eventPoint, out);
      m_found_intersection = true;
    }
  }

  
#ifndef NDEBUG
  void PrintEventQueue();
  void PrintSubCurves();
  void PrintStatusLine();
#endif

protected:

  /*! a pointer to a traits object */
  Traits *m_traits;

  Sweep_line_traits<Traits> m_sweep_line_traits;

  /*! an indication to whether the traits should be deleted in the distructor
   */
  bool m_traitsOwner;

  /*! if false, overlapping subcurves are reported only one. 
    Otherwise, they are reported as many times as they appeard. */
  bool m_overlapping;

  /*! if true, end points are reported as intersection points */
  bool m_includeEndPoints;

  /*! this is true when get_intersection_points() is called */
  bool m_stop_at_first_int;

  /*! the queue of events (intersection points) to handle */
  EventQueue *m_queue;

  /*! The subcurves, as created on the fly */
  Subcurve *m_subCurves;

  /*! The status line */
  StatusLine *m_statusLine;

  /*! a struct that holds params associated with the curve compare functor */
  CompareParams *m_comp_param;

  /*! The current position (in X) of the sweep line */
  Point_2 m_sweepLinePos;

  /*! if non x-monotone are specified, this hold the x-monotone 
    curves created when splitting them into x-monotone curves. */
  std::vector<X_monotone_curve_2> m_xcurves;

  /*! a pointer to the current event */
  Event *m_currentEvent;

  /*! when an intersection point is found this is turned to true */
  bool m_found_intersection;

  /*! An iterator of the  status line that is used as a hint for inserts. */
  StatusLineIter m_status_line_insert_hint;


  SubCurveList m_tmpOut;
  Point_2 m_lastReportedPoint;
  bool is_first_point;

  /*! An allocator for the events objects */
  EventAlloc m_eventAlloc;

  /*! An allocator for the SubCurve objects */
  SubCurveAlloc m_subCurveAlloc;

  /*! a master Event (created by default c'tor)to be used by allocator */
  Event m_masterEvent;

  /*! a master Subcurve (created by default c'tor) to be used by allocator */
  Subcurve m_masterSubcurve;


  /*! The num of subcurves  */
  unsigned int m_num_of_subCurves;

  

  
#ifndef NDEBUG
  int m_eventId;
#endif

  bool CurveStartsAtCurve(Subcurve *one, Subcurve *two)
  {
    SL_DEBUG(std::cout << "CurveStartsAtCurve: \n";)
    SL_DEBUG(std::cout << one->get_curve() << "\n" << two->get_curve() 
                       << "\n";)

    if ( m_traits->point_equal(one->get_left_end(),two->get_left_end()))
      return false;

    if ( !m_traits->point_equal(one->get_left_end(),
                                m_currentEvent->get_point()) )
      return false;
    if ( m_traits->curve_compare_y_at_x(one->get_left_end(), 
                                        two->get_curve()) == EQUAL )
      return true;
    return false;
  }

};

template <class CurveInputIterator,  class SweepLineTraits_2,
          class SweepEvent, class CurveWrap ,
          typename EventAlloc ,
          typename SubCurveAlloc >
Sweep_line_2_impl<CurveInputIterator,SweepLineTraits_2,SweepEvent,CurveWrap,
                  EventAlloc,SubCurveAlloc>::
~Sweep_line_2_impl() 
{
  if ( m_traitsOwner ) delete m_traits;
  
  for(unsigned int i=0 ; i < m_num_of_subCurves; ++i)
    m_subCurveAlloc.destroy(m_subCurves+i);

  if(m_num_of_subCurves) //if its zero, nothing to deallocate
    m_subCurveAlloc.deallocate(m_subCurves,m_num_of_subCurves); // deallocate memory 
  delete m_queue;
  delete m_statusLine;
  delete m_comp_param;
}



/*! initializes the data structures to work with:
 *  - x-monotonize the inf\put curves
 *  - for each end point of each curve create an event
 *  - initialize the event queue
 *  -
 */
template <class CurveInputIterator,  class SweepLineTraits_2,
          class SweepEvent, class CurveWrap,
          typename EventAlloc ,
          typename SubCurveAlloc >
inline void 
Sweep_line_2_impl<CurveInputIterator,SweepLineTraits_2,SweepEvent,CurveWrap,
                  EventAlloc,SubCurveAlloc>::
init(CurveInputIterator begin, CurveInputIterator end)
{
  PointLess pred(m_traits); // functor of the event queue(implemented as map)
  m_queue = new EventQueue(pred);  // allocate event queue (map)
  m_comp_param = new CompareParams(m_traits);
  StatusLineCurveLess slcurveless(m_comp_param , &m_sweepLinePos);//the functor of the status line (set)
  m_statusLine = new StatusLine(slcurveless); // allocate the status line
  m_status_line_insert_hint = m_statusLine->begin();

#ifndef NDEBUG
  m_eventId = 0;
#endif
  
  CurveInputIterator iter;
  for ( iter = begin ; iter != end ; ++iter)
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



/*! Given an x-monotone curve, create events for each end (if 
 *  one doesn't exist already). 
 *  For each curve create a Subcurve instance.
 */
template <class CurveInputIterator,  class SweepLineTraits_2,
          class SweepEvent, class CurveWrap,
          typename EventAlloc ,
          typename SubCurveAlloc >
inline void 
Sweep_line_2_impl<CurveInputIterator,SweepLineTraits_2,SweepEvent,CurveWrap,
                  EventAlloc,SubCurveAlloc>::
init_curve(X_monotone_curve_2 &curve,unsigned int j)
{ 
  Event *e ;
 
  m_subCurveAlloc.construct(m_subCurves+j,m_masterSubcurve);
  (m_subCurves+j)->init(curve);

  const Point_2 &left_end = (m_subCurves+j)->get_left_end();
  const Point_2 &right_end = (m_subCurves+j)->get_right_end();


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


    // handle the left end of the curve
    const std::pair<EventQueueIter, bool>& insertion_res2 =
      m_queue->insert(EventQueueValueType(left_end,0));  //the insertion return value
     
    if(insertion_res2.second == true)
    {
      e =  m_eventAlloc.allocate(1); // allocation of space for one dlink 
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
    e->add_curve_to_right(m_subCurves+j);
    PRINT_NEW_EVENT(left_end, e);
  
}




/*!
 * Perform intersection between the specified curve and all curves in the 
 * given group of curves.
 */ 
template <class CurveInputIterator, class SweepLineTraits_2,
          class SweepEvent, class CurveWrap,
          typename EventAlloc ,
          typename SubCurveAlloc >
inline void
Sweep_line_2_impl<CurveInputIterator,SweepLineTraits_2,SweepEvent,CurveWrap,
                  EventAlloc,SubCurveAlloc>::
intersect_curve_group(Subcurve *c1, SubCurveList &mylist)
{
  m_tmpOut.clear();
  SL_DEBUG(std::cout << "intersecting with " << mylist.size() << " curves\n";)
  SubCurveListIter i = mylist.begin();
  while ( i != mylist.end())
  {

    intersect(c1, *i);
    ++i;
  }
}



/*!
 * When a curve is removed from the status line for good, its top and
 * bottom neighbors become neighbors. This method finds these cases and
 * looks for the intersection point, if one exists.
 *
 * @param leftCurve a pointer to the curve that is about to be deleted
 * @return an iterator to the position where the curve will be removed from.
 */
template <class CurveInputIterator, class SweepLineTraits_2,
          class SweepEvent, class CurveWrap,
          typename EventAlloc ,
          typename SubCurveAlloc >
inline void
Sweep_line_2_impl<CurveInputIterator,SweepLineTraits_2,SweepEvent,CurveWrap,
                  EventAlloc,SubCurveAlloc>::
remove_curve_from_status_line(Subcurve *leftCurve)
{
  SL_DEBUG(std::cout << "remove_curve_from_status_line\n";)
  SL_DEBUG(PrintStatusLine();)
  SL_DEBUG(leftCurve->Print();)

  StatusLineIter sliter;
  sliter = leftCurve->get_hint(); 

  m_status_line_insert_hint = sliter; ++m_status_line_insert_hint;
  if ( !leftCurve->is_end_point(m_currentEvent->get_point())) {
    m_statusLine->erase(sliter);
    SL_DEBUG(std::cout << "remove_curve_from_status_line Done\n";)
    return;
  }

  CGAL_assertion(sliter!=m_statusLine->end());
  StatusLineIter end = m_statusLine->end(); --end;
  if ( sliter != m_statusLine->begin() && sliter != end ) 
  {
    SubCurveList mylist;
    StatusLineIter prev = sliter; --prev;
    
    // collect all curves that overlap with *prev
    StatusLineIter tmp = prev;
    mylist.push_back(*prev);
    while ( tmp != m_statusLine->begin() ) 
    {
      --tmp;
      if ( do_curves_overlap(*prev, *tmp))
        mylist.push_back(*tmp);
      else
        break;
    }
    
    StatusLineIter next = sliter; ++next;
    
    // intersect *next with the the *prev curve and all overlaps
    tmp = next;
    intersect_curve_group(*tmp, mylist);

    // if there are curves that overlap with the *next curve, intersect
    // them with the *prev curve and all overlaps
    ++tmp;
    while ( tmp != m_statusLine->end() ) 
    {
      if ( do_curves_overlap(*next, *tmp))
      {
        intersect_curve_group(*tmp, mylist);
        ++tmp;
      }
      else 
        break;
    }
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


template <class CurveInputIterator, class SweepLineTraits_2,
          class SweepEvent, class CurveWrap,
          typename EventAlloc ,
          typename SubCurveAlloc >
inline bool 
Sweep_line_2_impl<CurveInputIterator,SweepLineTraits_2,SweepEvent,CurveWrap,
                  EventAlloc,SubCurveAlloc>::
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
      e->add_curve_to_left(c2); 

      e->add_curve_to_right(c1);
      e->add_curve_to_right(c2);

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
              e->add_curve_to_right(c1);
        e->mark_internal_intersection_point();
      }

      e->add_curve_to_left(c2);
      if ( !c2->is_end_point(xp) ) 
      {
              e->add_curve_to_right(c2);
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




template <class CurveInputIterator,  class SweepLineTraits_2,
          class SweepEvent, class CurveWrap,
          typename EventAlloc ,
          typename SubCurveAlloc >
inline bool
Sweep_line_2_impl<CurveInputIterator,SweepLineTraits_2,SweepEvent,CurveWrap,
                  EventAlloc,SubCurveAlloc>::
do_curves_overlap(Subcurve *c1, Subcurve *c2)
{
  SL_DEBUG(std::cout << "do_curves_overlap " << m_sweepLinePos << "\n" 
                     << "\t" << c1->get_curve() << "\n"
                     << "\t" << c2->get_curve() << "\n";)

  const Point_2 *p = &(c2->get_last_point());
  if ( m_traits->compare_xy(c1->get_last_point(), 
                           c2->get_last_point()) == LARGER ) 
    p = &(c1->get_last_point());

  if ((m_traits->curves_compare_y_at_x(c1->get_curve(),
                                       c2->get_curve(),
                                       *p) != EQUAL))
    return false;

  if ( m_traits->curves_overlap(c1->get_curve(),c2->get_curve()) )
    return true;

  return false;
}

template <class CurveInputIterator,  class SweepLineTraits_2,
          class SweepEvent, class CurveWrap,
          typename EventAlloc,
          typename SubCurveAlloc >
inline bool
Sweep_line_2_impl<CurveInputIterator,SweepLineTraits_2,SweepEvent,CurveWrap,
                  EventAlloc,SubCurveAlloc>::
similar_curves(const X_monotone_curve_2 &a, const X_monotone_curve_2 &b)
{
  return ( m_traits->curve_equal(a, b));
}




////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//                         DEBUG UTILITIES                                //
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

#ifndef NDEBUG

template <class CurveInputIterator,  class SweepLineTraits_2,
          class SweepEvent, class CurveWrap,
          typename EventAlloc ,
          typename SubCurveAlloc>
inline void 
Sweep_line_2_impl<CurveInputIterator,SweepLineTraits_2,SweepEvent,CurveWrap,
                  EventAlloc,SubCurveAlloc>::
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

template <class CurveInputIterator,  class SweepLineTraits_2,
          class SweepEvent, class CurveWrap,
          typename EventAlloc ,
          typename SubCurveAlloc >
inline void 
Sweep_line_2_impl<CurveInputIterator,SweepLineTraits_2,SweepEvent,CurveWrap,
                  EventAlloc,SubCurveAlloc>::
PrintSubCurves()
{
  SL_DEBUG(std::cout << std::endl << "Sub curves: " << std::endl;)
  for(unsigned int i=0 ; i < m_num_of_subCurves ; ++i)
  {
    m_subCurves[i].Print();
  }
}

template <class CurveInputIterator,  class SweepLineTraits_2,
          class SweepEvent, class CurveWrap,
          typename EventAlloc ,
          typename SubCurveAlloc >
inline void 
Sweep_line_2_impl<CurveInputIterator,SweepLineTraits_2,SweepEvent,CurveWrap,
                  EventAlloc,SubCurveAlloc>::
PrintStatusLine()
{
  if ( m_statusLine->size() == 0) {
    std::cout << std::endl << "Status line: empty" << std::endl;
    return;
  }
  std::cout << std::endl << "Status line: (" 
            << m_sweepLinePos << ")" << std::endl;
  StatusLineIter iter = m_statusLine->begin();
  while ( iter != m_statusLine->end() )
  {
    (*iter)->Print();
    ++iter;
  }
  std::cout << "Status line - end" << std::endl;
}


#endif

CGAL_END_NAMESPACE

#endif
