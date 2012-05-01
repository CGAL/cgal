// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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

/*! \file
 * Definition of the Basic_sweep_line_2 class.
 */

#include <boost/mpl/assert.hpp>
#include <CGAL/assertions.h>
#include <CGAL/memory.h>
#include <CGAL/Sweep_line_2/Sweep_line_functors.h>
#include <CGAL/Sweep_line_2/Sweep_line_subcurve.h>
#include <CGAL/Sweep_line_2/Sweep_line_event.h>
#include <CGAL/Multiset.h>
#include <CGAL/Arrangement_2/Arr_traits_adaptor_2.h>
#include <CGAL/Arr_tags.h>
#include <vector>
#include <algorithm>
#include <iterator>

#ifndef CGAL_SL_VERBOSE

#define CGAL_SL_DEBUG(a)
#define CGAL_PRINT_INSERT(a)
#define CGAL_PRINT_ERASE(a)
#define CGAL_PRINT_NEW_EVENT(p, e) 
#define CGAL_PRINT_UPDATE_EVENT(p, e) 
#define CGAL_PRINT(a)

#else

#include <iostream>

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
#define CGAL_PRINT_UPDATE_EVENT(p, e) \
{ std::cout << "%%% an event was updated at " << (p) << std::endl; \
  (e)->Print(); }
#define CGAL_PRINT(a) { std::cout << a ; }

#endif



namespace CGAL {

/*! \class Basic_Sweep_line_2 
 * A class that implements the sweep line algorithm for general x-monotone
 * curves that are pairwise disjoint in their interiors (an additional set
 * of isolated points may also be supplied).
 * The x-montone curve type and the point type are defined by the traits class
 * that is one of the template parameters.
 */
template < class Traits_,
           class Visitor_,
           class Subcurve_ = Sweep_line_subcurve<Traits_>,
           typename Event_ = Sweep_line_event<Traits_, Subcurve_>,
           typename Allocator_ = CGAL_ALLOCATOR(int) >
class Basic_sweep_line_2
{
public:

  typedef Traits_                                       Traits_2;
  typedef Visitor_                                      Visitor;
  typedef Event_                                        Event;
  typedef Subcurve_                                     Subcurve;
  typedef Allocator_                                    Allocator;

  typedef Arr_traits_basic_adaptor_2<Traits_2>          Traits_adaptor_2;
  typedef typename Traits_adaptor_2::Point_2            Point_2;
  typedef typename Traits_adaptor_2::X_monotone_curve_2 X_monotone_curve_2;

  typedef typename Traits_adaptor_2::Left_side_category   Left_side_category;
  typedef typename Traits_adaptor_2::Bottom_side_category Bottom_side_category;
  typedef typename Traits_adaptor_2::Top_side_category    Top_side_category;
  typedef typename Traits_adaptor_2::Right_side_category  Right_side_category;

  BOOST_MPL_ASSERT(
      (typename 
       Arr_sane_identified_tagging< Left_side_category, Bottom_side_category, 
       Top_side_category, Right_side_category >::result)
  );
  
protected:

  typedef typename Arr_are_all_sides_oblivious_tag< 
                     Left_side_category, Bottom_side_category, 
                     Top_side_category, Right_side_category >::result
  Are_all_sides_oblivious_tag;
  
public:

  typedef CGAL::Compare_events<Traits_adaptor_2, Event> Compare_events;
  typedef Multiset<Event*, Compare_events, Allocator>   Event_queue; 
  typedef typename Event_queue::iterator                Event_queue_iterator;

  typedef typename Event::Subcurve_iterator
    Event_subcurve_iterator;

  typedef Sweep_line_event<Traits_2, Subcurve>          Base_event;
  typedef typename Base_event::Attribute                Attribute;
  
  typedef Sweep_line_subcurve<Traits_2>                 Base_subcurve;
  typedef class Curve_comparer<Traits_2, Base_subcurve> Compare_curves;
  typedef Multiset<Base_subcurve*,
                   Compare_curves, 
                   Allocator>                           Status_line;
  typedef typename Status_line::iterator                Status_line_iterator;

  typedef typename Allocator::template rebind<Event>    Event_alloc_rebind;
  typedef typename Event_alloc_rebind::other            Event_alloc;

  typedef typename Allocator::template rebind<Subcurve> Subcurve_alloc_rebind;
  typedef typename Subcurve_alloc_rebind::other         Subcurve_alloc;


protected:

  /*! \struct
   * An auxiliary functor for comparing event pointers.
   */
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

  // Data members:
  const Traits_adaptor_2 * m_traits;// A traits-class object.
  bool              m_traitsOwner;  // Whether this object was allocated by
                                    // this class (and thus should be freed).

  Event            *m_currentEvent; // The current event.

  Compare_curves   m_statusLineCurveLess;
                                     // Comparison functor for the status line.

  Compare_events   m_queueEventLess; // Comparison functor for the event queue.

  Event_queue     *m_queue;          // The event queue (the X-structure).

  Subcurve        *m_subCurves;      // An array of the subcurves.
  Status_line      m_statusLine;     // The status line (the Y-structure).

  Allocated_events_set m_allocated_events;
                                     // The events that have been allocated
                                     // (and have not yet been deallocated).

  Status_line_iterator m_status_line_insert_hint;
                                     // An iterator of the status line, which
                                     // is used as a hint for insertions.

  bool             m_is_event_on_above;
                                     // Indicates if the current event is on
                                     // the interior of existing curve. This 
                                     // may happen only with events that are
                                     // associated with isolated query points.

  Event_alloc    m_eventAlloc;       // An allocator for the events objects.
  Subcurve_alloc m_subCurveAlloc;    // An allocator for the subcurve objects.

  Event          m_masterEvent;      // A master Event (created once by the
                                     // constructor) for the allocator's usage.

  Subcurve       m_masterSubcurve;   // A master Subcurve (created once by the
                                     // constructor) for the allocator's usage.

  unsigned int   m_num_of_subCurves; // Number of subcurves.

  Visitor       *m_visitor;          // The sweep-line visitor that will be
                                     // notified during the sweep.

public:

  /*!
   * Constructor.
   * \param visitor A pointer to a sweep-line visitor object.
   */
  Basic_sweep_line_2 (Visitor *visitor);

  /*!
   * Constructor with a traits class.
   * \param traits A pointer to a sweep-line traits object.
   * \param visitor A pointer to a sweep-line visitor object.
   */
  Basic_sweep_line_2 (const Traits_2 *traits, Visitor *visitor);

  /*! Destrcutor. */
  virtual ~Basic_sweep_line_2 ();

  /*!
   * Run the sweep-line algorithm on a given range of x-monotone curves.
   * \param curves_begin An iterator for the first curve in the range.
   * \param curves_end A past-the-end iterator for the range.
   * \pre The value-type of CurveInputIterator is X_monotone_curve_2.
   */
  template<class CurveInputIterator>
  void sweep (CurveInputIterator curves_begin,
              CurveInputIterator curves_end)
  {
    m_visitor->before_sweep();
    _init_sweep(curves_begin, curves_end);
    //m_visitor ->after_init();
    _sweep();
    _complete_sweep();
    m_visitor ->after_sweep();
  }

  /*!
   * Run the sweep-line algorithm on a range of x-monotone curves and a range 
   * of action event points (if a curve passed through an action point, it will
   * be split).
   * \param curves_begin  An iterator for the first x-monotone curve in the
   *                      range.
   * \param curves_end A past-the-end iterator for this range.
   * \param points_begin An iterator for the first point in the range.
   * \param points_end A past-the-end iterator for this range.
   * \pre The value-type of XCurveInputIterator is the traits-class
   *      X_monotone_curve_2, and the value-type of PointInputIterator is the
   *      traits-class Point_2.
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
    //m_visitor ->after_init();
    _sweep();
    _complete_sweep();
    m_visitor ->after_sweep();
  }

  /*!
   * Run the sweep-line alogrithm on a range of x-monotone curves, a range   
   * of action event points (if a curve passed through an action point, it will
   * be split) and a range of query points (if a curve passed through a
   * query point,it will not be splitted).
   * \param curves_begin An iterator for the first x-monotone curve in the
   *                     range.
   * \param curves_end A past-the-end iterator for this range.
   * \param points_begin An iterator for the first point in the range.
   * \param points_end A past-the-end iterator for this range.
   * \pre The value-type of XCurveInputIterator is the traits-class 
   *      X_monotone_curve_2, and the value-type of PointInputIterator is the 
   *      traits-class Point_2.
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
    //m_visitor ->after_init();
    _sweep();
    _complete_sweep();
    m_visitor ->after_sweep();
  }

  /*! Get an iterator for the first subcurve in the status line. */
  Status_line_iterator status_line_begin ()
  {
    return (m_statusLine.begin());
  }

  /*! Get a past-the-end iterator for the subcurves in the status line. */
  Status_line_iterator status_line_end()
  {
    return (m_statusLine.end());
  }

  /*! Get the status line size. */
  unsigned int status_line_size() const
  {
    return (m_statusLine.size());
  }

  /*! Check if the status line is empty. */
  bool is_status_line_empty() const
  {
    return (m_statusLine.empty());
  }

  /*! Get an iterator for the first event in event queue. */
  Event_queue_iterator event_queue_begin()
  {
    return (m_queue->begin());
  }

  /*! Get a past-the-end iterator for the events in the in event queue. */
  Event_queue_iterator event_queue_end()
  {
    return (m_queue->end());
  }

   /*! Get the event queue size. */
  unsigned int event_queue_size() const
  {
    return (m_queue->size());
  }

  /*! Check if the event queue is empty. */
  bool is_event_queue_empty() const
  {
    return (m_queue->empty());
  }

  /*! 
   * Stop the sweep by erasing the event queue (except for the current event).
   * This function may called by the visitor during 'arter_handle_event' in
   * order to stop the sweep-line process.
   */
  void stop_sweep();

  /*!
   * Deallocate event object.
   * This method is made public to allow the visitor to manage the events
   * deallocation (as necessary). 
   */
  void deallocate_event(Event* event);

  /*! Get the current event */
  Event* current_event()
  {
    return (m_currentEvent);
  }

  /*! Get the traits object */
  const Traits_2 * traits ()
  {
    return m_traits;
  }

protected:

  /*! Perform the main sweep-line loop. */
  void _sweep();

  /*! Create an event object for each input point. */
  template <class PointInputIterator>
  void _init_points (PointInputIterator points_begin,
                     PointInputIterator points_end,
                     Attribute type)
  {
    PointInputIterator   pit;
    for (pit = points_begin; pit != points_end; ++pit)
      _init_point (*pit, type);

    return;
  }

  /*! Create a Subcurve object and two Event objects for each curve. */
  template<class CurveInputIterator>
  void _init_curves (CurveInputIterator curves_begin,
                     CurveInputIterator curves_end)
  {
    CurveInputIterator   cit;
    unsigned int         index = 0;

    for (cit = curves_begin; cit != curves_end; ++cit, ++index)
      _init_curve (*cit, index);

    return;
  }

  /*! Initiliaze the sweep algorithm. */
  template<class CurveInputIterator>
  void _init_sweep (CurveInputIterator curves_begin,
                    CurveInputIterator curves_end)
  {
    // m_num_of_subCurves should be a size_t for "huge" data sets
    m_num_of_subCurves = static_cast<int>(std::distance (curves_begin, curves_end));

    _init_structures();

    // Initialize the curves.
    _init_curves (curves_begin, curves_end);
    return;
  }

  /*! Initialize the data structures for the sweep-line algorithm. */
  virtual void _init_structures ();

  /*! Compete the sweep (compete data strcures) */
  virtual void _complete_sweep();

  /*!
   * Initialize an event associated with a point.
   * \param p The given point.
   * \param type The event type.
   */
  void _init_point (const Point_2& pt, Attribute type);

  /*!
   * Initialize the events associated with an x-monotone curve.
   * \param curve The given x-monotone curve.
   * \param index Its unique index.
   */
  void _init_curve (const X_monotone_curve_2& curve, unsigned int index);

  /*!
   * Initialize an event associated with an x-monotone curve end.
   * \param cv The given x-monotone curve.
   * \param ind Its end (ARR_MIN_END or ARR_MAX_END).
   * \param sc The subcurve corresponding to cv.
   */
  void _init_curve_end (const X_monotone_curve_2& cv,
                        Arr_curve_end ind,
                        Subcurve* sc);
  
  /*!
   * Handle the subcurves that are to the left of the event point (i.e., 
   * subcurves that we are done with).
   */
  virtual void _handle_left_curves();

  /*!
   * Handle an event that does not have any incident left curves.
   * Such an event is usually the left endpoint of its incident right
   * subcurves, and we locate thei position in the status line.
   */
  void _handle_event_without_left_curves ();

  /*!
   * Sort the left subcurves of an event point according to their order in
   * their status line (no geometric comprasions are needed).
   */
  void _sort_left_curves ();

  /*! Handle the subcurves to the right of the current event point. */
  virtual void _handle_right_curves ();

  /*!
   * Add a subcurve to the right of an event point.
   * \param event The event point.
   * \param curve The subcurve to add.
   * \return (true) if an overlap occured; (false) otherwise.
   */
  virtual bool _add_curve_to_right (Event* event, Subcurve* curve,
                                    bool overlap_exist = false);

  /*! Remove a curve from the status line. */
  void _remove_curve_from_status_line (Subcurve *leftCurve);
 
  /*!
   * Allocate an event object associated with a given point.
   * \param pt The point.
   * \param type The event type.
   * \param ps_x The location of the point in x.
   * \param ps_y The location of the point in y.
   * \pre Neither one of the boundary conditions is +/-oo. 
   * \return The created event.
   */
  Event* _allocate_event (const Point_2& pt, Attribute type,
                          Arr_parameter_space ps_x, Arr_parameter_space ps_y);

  /*!
   * Allocate an event at open boundary, 
   * which is not associated with a valid point.
   * \param type The event type.
   * \param ps_x The location of the point in x.
   * \param ps_y The location of the point in y.
   * \param At least one of the boundary conditions is +/-oo.
   * \return The created event.
   */
  Event* _allocate_event_at_open_boundary (Attribute type,
                                           Arr_parameter_space ps_x,
                                           Arr_parameter_space ps_y);

  /*! 
   * Push a finite event point into the event queue.
   * \param pt The point associated with the event.
   * \param type The event type.
   * \param ps_x The location of the point in x.
   * \param ps_y The location of the point in y.
   * \param sc A subcurve that the new event represents on of its endpoints.
   * \return A pair that comprises a pointer to the event, and a flag
   *         indicating whether this is a new event (if false, the event
   *         was in the queue and we just updated it).
   */
  std::pair<Event*, bool> _push_event (const Point_2& pt,
                                       Attribute type,
                                       Arr_parameter_space ps_x,
                                       Arr_parameter_space ps_y,
                                       Subcurve* sc = NULL);

  /*! 
   * Push an event point associated with a curve end into the event queue.
   * \param cv The x-monotone curve.
   * \param ind The relevant curve end.
   * \param type The event type.
   * \param ps_x The location of the point in x.
   * \param ps_y The location of the point in y.
   * \param sc A subcurve that the new event represents on of its endpoints.
   * \return A pair that comprises a pointer to the event, and a flag
   *         indicating whether this is a new event (if false, the event
   *         was in the queue and we just updated it).
   */
  std::pair<Event*, bool> _push_event (const X_monotone_curve_2& cv,
                                       Arr_curve_end ind,
                                       Attribute type,
                                       Arr_parameter_space ps_x,
                                       Arr_parameter_space ps_y,
                                       Subcurve* sc = NULL);

  void _update_event_at_open_boundary(Event* e,
                                      const X_monotone_curve_2& cv,
                                      Arr_curve_end ind,
                                      bool is_new)
  {
    _update_event_at_open_boundary(e, cv, ind, is_new, 
                                   Are_all_sides_oblivious_tag());
  }

  void _update_event_at_open_boundary(Event* e,
                                      const X_monotone_curve_2& cv,
                                      Arr_curve_end ind,
                                      bool is_new,
                                      Arr_not_all_sides_oblivious_tag)
  {
    m_visitor->update_event (e, cv, ind, is_new);
  }

  void _update_event_at_open_boundary(Event* /* e */,
                                      const X_monotone_curve_2& /* cv */,
                                      Arr_curve_end /* ind */,
                                      bool /* is_new */,
                                      Arr_all_sides_oblivious_tag)
  {
    CGAL_error();
  }

#ifdef CGAL_SL_VERBOSE
  void PrintEventQueue();
  void PrintSubCurves();
  void PrintStatusLine();
  void PrintOpenBoundaryType(Arr_parameter_space x, Arr_parameter_space y);
  void PrintEvent(const Event* e);
#endif

};

//DEBUG UTILITIES
#ifdef CGAL_SL_VERBOSE
  #include <CGAL/Sweep_line_2/Sweep_line_2_debug.h>
#endif

} //namespace CGAL

#include <CGAL/Sweep_line_2/Basic_sweep_line_2_impl.h>

#endif
