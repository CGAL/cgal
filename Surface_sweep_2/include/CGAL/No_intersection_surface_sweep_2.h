// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Baruch Zukerman <baruchzu@post.tau.ac.il>
//             Efi Fogel <efifogel@gmail.com>
//               (based on old version by Tali Zvi)

#ifndef CGAL_NO_INTERSECTION_SURFACE_SWEEP_2_H
#define CGAL_NO_INTERSECTION_SURFACE_SWEEP_2_H

/*! \file
 *
 * Definition of the No_intersection_surface_sweep_2 class.
 */

#include <boost/mpl/assert.hpp>

#include <CGAL/license/Surface_sweep_2.h>
#include <CGAL/assertions.h>
#include <CGAL/memory.h>
#include <CGAL/Surface_sweep_2/Event_comparer.h>
#include <CGAL/Surface_sweep_2/Curve_comparer.h>
#include <CGAL/Multiset.h>
#include <CGAL/Arrangement_2/Arr_traits_adaptor_2.h>
#include <CGAL/Arr_tags.h>
#include <vector>
#include <algorithm>
#include <iterator>

#ifndef CGAL_SS_VERBOSE

#define CGAL_SS_DEBUG(a)
#define CGAL_SS_PRINT(a)
#define CGAL_SS_PRINT_TEXT(text)
#define CGAL_SS_PRINT_EOL()
#define CGAL_SS_PRINT_START(name)
#define CGAL_SS_PRINT_START_EOL(name)
#define CGAL_SS_PRINT_END(name)
#define CGAL_SS_PRINT_END_EOL(name)
#define CGAL_SS_PRINT_CURVE(text)
#define CGAL_SS_PRINT_EVENT_INFO(event)
#define CGAL_SS_PRINT_STATUS_LINE()

#define CGAL_SS_PRINT_INSERT(a)
#define CGAL_SS_PRINT_ERASE(a)
#define CGAL_SS_PRINT_NEW_EVENT(p, e)
#define CGAL_SS_PRINT_UPDATE_EVENT(p, e)

#else

#include <iostream>

#define CGAL_SS_DEBUG(a) {a;}
#define CGAL_SS_PRINT(a) std::cout << a
#define CGAL_SS_PRINT_TEXT(text) this->print_text(text)
#define CGAL_SS_PRINT_EOL() this->print_eol()
#define CGAL_SS_PRINT_START(name) this->print_start(name, false)
#define CGAL_SS_PRINT_START_EOL(name) this->print_start(name, true)
#define CGAL_SS_PRINT_END(name) this->print_end(name, false)
#define CGAL_SS_PRINT_END_EOL(name) this->print_end(name, true)
#define CGAL_SS_PRINT_CURVE(text) this->print_curve(text)
#define CGAL_SS_PRINT_EVENT_INFO(event) this->print_event_info(event)
#define CGAL_SS_PRINT_STATUS_LINE() this->PrintStatusLine()

#define CGAL_SS_PRINT_INSERT(a) {           \
    this->print_text("+++ inserting ");     \
    (a)->Print();                           \
    this->print_eol();                      \
    this->print_text(" currentPos = ");     \
    this->PrintEvent(this->m_currentEvent); \
    this->print_eol();                      \
  }
#define CGAL_SS_PRINT_ERASE(a) {      \
    this->print_text("--- erasing "); \
    (a)->Print();                     \
    this->print_eol();                \
  }
#define CGAL_SS_PRINT_NEW_EVENT(p, e) {                  \
    this->print_text("%%% a new event was created at "); \
    CGAL_SS_PRINT(p);                                    \
    this->print_eol();                                   \
  }
#define CGAL_SS_PRINT_UPDATE_EVENT(p, e) {            \
    this->print_text("%%% an event was updated at "); \
    CGAL_SS_PRINT(p);                                 \
    this->print_eol();                                \
    this->print_event_info(e);                        \
  }

#endif

namespace CGAL {
namespace Surface_sweep_2 {

/*! \class No_intersection_surface_sweep_2
 * A class that implements the sweep line algorithm for general x-monotone
 * curves that are pairwise disjoint in their interiors (an additional set
 * of isolated points may also be supplied).
 * The x-montone curve type and the point type are defined by the traits class
 * that is one of the template parameters.
 */
template <typename Visitor_>
class No_intersection_surface_sweep_2 {
public:
  typedef Visitor_                                      Visitor;

  typedef typename Visitor::Geometry_traits_2           Geometry_traits_2;
  typedef typename Visitor::Event                       Event;
  typedef typename Visitor::Subcurve                    Subcurve;
  typedef typename Visitor::Allocator                   Allocator;

private:
  typedef Geometry_traits_2                             Gt2;

public:
  typedef Arr_traits_basic_adaptor_2<Gt2>               Traits_adaptor_2;
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
    Are_all_sides_oblivious_category;

public:
  typedef CGAL::Surface_sweep_2::Event_comparer<Traits_adaptor_2, Event>
                                                        Event_comparer;
  typedef Multiset<Event*, Event_comparer, Allocator, Tag_true>
                                                        Event_queue;
  typedef typename Event_queue::iterator                Event_queue_iterator;

  typedef typename Event::Subcurve_iterator
    Event_subcurve_iterator;
  typedef typename Event::Subcurve_const_iterator
    Event_subcurve_const_iterator;

  typedef typename Event::Attribute                     Attribute;

  typedef class Curve_comparer<Gt2, Event, Subcurve>    Compare_curves;
  typedef Multiset<Subcurve*, Compare_curves, Allocator>
                                                        Status_line;
  typedef typename Status_line::iterator                Status_line_iterator;

  typedef std::allocator_traits<Allocator> Allocator_traits;
  typedef typename Allocator_traits::template rebind_alloc<Subcurve> Subcurve_alloc;

protected:

  typedef Compact_container<Event>                 Allocated_events_set;
  typedef typename Allocated_events_set::iterator  Allocated_events_iterator;

  // Data members:
  const Traits_adaptor_2* m_traits; // A traits-class object.
  bool m_traitsOwner;               // Whether this object was allocated by
                                    // this class (and thus should be freed).

  Event* m_currentEvent;            // The current event.

  Compare_curves m_statusLineCurveLess;
                                    // Comparison functor for the status line.

  Event_comparer m_queueEventLess;  // Comparison functor for the event queue.

  Event_queue* m_queue;             // The event queue (the X-structure).

  Subcurve* m_subCurves;            // An array of the subcurves.
  Status_line m_statusLine;         // The status line (the Y-structure).

  Allocated_events_set m_allocated_events;
                                    // The events that have been allocated
                                    // (and have not yet been deallocated).

  Status_line_iterator m_status_line_insert_hint;
                                    // An iterator of the status line, which
                                    // is used as a hint for insertions.

  bool m_is_event_on_above;         // Indicates if the current event is on
                                    // the interior of existing curve. This
                                    // may happen only with events that are
                                    // associated with isolated query points.

  Subcurve_alloc m_subCurveAlloc;   // An allocator for the subcurve objects.

  Event m_masterEvent;              // A master Event (created once by the
                                    // constructor) for the allocator's usage.

  Subcurve m_masterSubcurve;        // A master Subcurve (created once by the
                                    // constructor) for the allocator's usage.

  //! \todo m_num_of_subCurves should be a size_t for "huge" data sets
  unsigned int m_num_of_subCurves;  // Number of subcurves.

  Visitor* m_visitor;               // The sweep-line visitor that will be
                                    // notified during the sweep.

public:
  /*! Constructor.
   * \param visitor A pointer to a sweep-line visitor object.
   */
  No_intersection_surface_sweep_2(Visitor* visitor);

  /*! Constructor with a traits class.
   * \param traits A pointer to a sweep-line traits object.
   * \param visitor A pointer to a sweep-line visitor object.
   */
  No_intersection_surface_sweep_2(const Gt2* traits, Visitor* visitor);

  /*! Destructor. */
  virtual ~No_intersection_surface_sweep_2();

  /*! Run the sweep-line algorithm on a given range of x-monotone curves.
   * \param curves_begin An iterator for the first curve in the range.
   * \param curves_end A past-the-end iterator for the range.
   * \pre The value-type of CurveInputIterator is X_monotone_curve_2.
   */
  template <typename CurveInputIterator>
  void sweep(CurveInputIterator curves_begin, CurveInputIterator curves_end)
  {
    m_visitor->before_sweep();
    _init_sweep(curves_begin, curves_end);
    //m_visitor->after_init();
    _sweep();
    _complete_sweep();
    m_visitor->after_sweep();
  }

  /*! Run the sweep-line algorithm on a range of x-monotone curves and a range
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
  template <typename CurveInputIterator, class PointInputIterator>
  void sweep(CurveInputIterator curves_begin,
             CurveInputIterator curves_end,
             PointInputIterator action_points_begin,
             PointInputIterator action_points_end)
  {
    m_visitor->before_sweep();
    _init_sweep(curves_begin, curves_end);
    _init_points(action_points_begin, action_points_end, Event::ACTION);
    //m_visitor->after_init();
    _sweep();
    _complete_sweep();
    m_visitor->after_sweep();
  }

  /*! Run the sweep-line alogrithm on a range of x-monotone curves, a range
   * of action event points (if a curve passed through an action point, it will
   * be split) and a range of query points (if a curve passed through a
   * query point,it will not be split).
   * \param curves_begin An iterator for the first x-monotone curve in the
   *                     range.
   * \param curves_end A past-the-end iterator for this range.
   * \param points_begin An iterator for the first point in the range.
   * \param points_end A past-the-end iterator for this range.
   * \pre The value-type of XCurveInputIterator is the traits-class
   *      X_monotone_curve_2, and the value-type of PointInputIterator is the
   *      traits-class Point_2.
   */
  template <typename CurveInputIterator, typename ActionPointItr,
            typename QueryPointItr>
  void sweep(CurveInputIterator curves_begin,
             CurveInputIterator curves_end,
             ActionPointItr action_points_begin,
             ActionPointItr action_points_end,
             QueryPointItr query_points_begin,
             QueryPointItr query_points_end)
  {
    m_visitor->before_sweep();
    _init_sweep(curves_begin, curves_end);
    _init_points(action_points_begin, action_points_end, Event::ACTION);
    _init_points(query_points_begin, query_points_end, Event::QUERY);
    //m_visitor->after_init();
    _sweep();
    _complete_sweep();
    m_visitor->after_sweep();
  }

  /*!
    runs the indexed sweep-line algorithm on a given range of
    x-monotone curves and accessor. The main difference from the
    original `sweep()` function is that the accessor allows to get
    indices of end-points and avoid checking several times in the
    event queue when the same vertex has several incident edges.

    \param edges A range of edges
    \param accessor An object providing, for a value_type `e` of
           range `edges`, the following methods:
     - size_t             min_end_index(e) -> the index of the min end of curve c
     - size_t             max_end_index(e) -> the index of the max end of curve c
     - X_monotone_curve_2 curve(e)         -> the x-monotone curve associated to e
     - size_t             nb_vertices()    -> the total number of points/events
     - void               before_init()    -> called before initialization
     - void               after_init()     -> called after initialization
   */
  template <typename EdgeRange, typename Accessor>
  void indexed_sweep (const EdgeRange& edges,
                      const Accessor& accessor)
  {
    m_visitor->before_sweep();
    accessor.before_init();
    _init_indexed_sweep(edges, accessor);
    accessor.after_init();
    _sweep();
    _complete_sweep();
    m_visitor->after_sweep();
  }

  /*!
    runs the indexed sweep-line algorithm on a given range of
    x-monotone curves and accessor. The main difference from the
    original `sweep()` function is that the accessor allows to get
    indices of end-points and avoid checking several times in the
    event queue when the same vertex has several incident edges.

    Variant with action event points (if a curve passed through an
    action point, it will be split).

    \param edges A range of edges
    \param accessor An object providing, for a value_type `e` of
           range `edges`, the following methods:
     - size_t             min_end_index(e) -> the index of the min end of curve c
     - size_t             max_end_index(e) -> the index of the max end of curve c
     - X_monotone_curve_2 curve(e)         -> the x-monotone curve associated to e
     - size_t             nb_vertices()    -> the total number of points/events
     - void               before_init()    -> called before initialization
     - void               after_init()     -> called after initialization
    \param points_begin An iterator for the first point in the range.
    \param points_end A past-the-end iterator for this range.
    \pre The value-type of PointInputIterator is the
         traits-class Point_2.
  */
  template <typename EdgeRange, typename Accessor,
            typename PointInputIterator>
  void indexed_sweep (const EdgeRange& edges,
                      const Accessor& accessor,
                      PointInputIterator action_points_begin,
                      PointInputIterator action_points_end)
  {
    m_visitor->before_sweep();
    accessor.before_init();
    _init_indexed_sweep(edges, accessor);
    accessor.after_init();
    _init_points(action_points_begin, action_points_end, Event::ACTION);
    _sweep();
    _complete_sweep();
    m_visitor->after_sweep();
  }

  /*! Get an iterator for the first subcurve in the status line. */
  Status_line_iterator status_line_begin() { return m_statusLine.begin(); }

  /*! Get a past-the-end iterator for the subcurves in the status line. */
  Status_line_iterator status_line_end() { return m_statusLine.end(); }

  /*! Get the status line size. */
  unsigned int status_line_size() const { return m_statusLine.size(); }

  /*! Check if the status line is empty. */
  bool is_status_line_empty() const { return m_statusLine.empty(); }

  /*! Get an iterator for the first event in event queue. */
  Event_queue_iterator event_queue_begin() { return m_queue->begin(); }

  /*! Get a past-the-end iterator for the events in the in event queue. */
  Event_queue_iterator event_queue_end() { return m_queue->end(); }

   /*! Get the event queue size. */
  unsigned int event_queue_size() const { return m_queue->size(); }

  /*! Check if the event queue is empty. */
  bool is_event_queue_empty() const { return m_queue->empty(); }

  /*! Stop the sweep by erasing the event queue (except for the current event).
   * This function may called by the visitor during 'arter_handle_event' in
   * order to stop the sweep-line process.
   */
  void stop_sweep();

  /*! Deallocate event object.
   * This method is made public to allow the visitor to manage the events
   * deallocation (as necessary).
   */
  void deallocate_event(Event* event);

  /*! Get the current event */
  Event* current_event() { return m_currentEvent; }

  /*! Get the traits object */
  const Gt2* traits() { return m_traits; }

protected:
  /*! Perform the main sweep-line loop. */
  void _sweep();

  /*! Create an event object for each input point. */
  template <typename PointInputIterator>
  void _init_points(PointInputIterator points_begin,
                    PointInputIterator points_end,
                    Attribute type)
  {
    for (PointInputIterator pit = points_begin; pit != points_end; ++pit)
      _init_point(*pit, type);
  }

  /*! Create a Subcurve object and two Event objects for each curve. */
  template <typename CurveInputIterator>
  void _init_curves(CurveInputIterator curves_begin,
                    CurveInputIterator curves_end)
  {
    CurveInputIterator cit;
    unsigned int index = 0;
    for (cit = curves_begin; cit != curves_end; ++cit, ++index)
      _init_curve(*cit, index);
  }

  /*! Create a Subcurve object and two Event objects for each curve. */
  template <typename EdgeRange, typename Accessor>
  void _init_indexed_curves(const EdgeRange& edges,
                            const Accessor& accessor)
  {
    std::vector<Event_queue_iterator> events (accessor.nb_vertices());

    unsigned int index = 0;
    for (const auto& e : edges)
    {
      std::size_t max_end = accessor.max_end_index(e);
      std::size_t min_end = accessor.min_end_index(e);
      const X_monotone_curve_2& curve = accessor.curve (e);

      // Construct and initialize a subcurve object.
      std::allocator_traits<Subcurve_alloc>::construct(m_subCurveAlloc, m_subCurves + index, m_masterSubcurve );
      (m_subCurves + index)->set_hint(this->m_statusLine.end());
      (m_subCurves + index)->init (curve);

      _init_curve_end(curve, ARR_MAX_END, m_subCurves + index, events, max_end);
      _init_curve_end(curve, ARR_MIN_END, m_subCurves + index, events, min_end);

      ++ index;
    }
  }

  /*! Initiliaze the sweep algorithm. */
  template <typename CurveInputIterator>
  void _init_sweep(CurveInputIterator curves_begin,
                   CurveInputIterator curves_end)
  {
    m_num_of_subCurves =
      static_cast<unsigned int>(std::distance(curves_begin, curves_end));
    _init_structures();
    _init_curves(curves_begin, curves_end);     // initialize the curves
  }

  /*! Initiliaze the sweep algorithm. */
  template <typename EdgeRange, typename Accessor>
  void _init_indexed_sweep(const EdgeRange& edges,
                           const Accessor& accessor)
  {
    m_num_of_subCurves =
      static_cast<unsigned int>(std::distance(edges.begin(), edges.end()));
    _init_structures();
    _init_indexed_curves(edges, accessor);     // initialize the curves
  }

  /*! Initialize the data structures for the sweep-line algorithm. */
  virtual void _init_structures();

  /*! Complete the sweep (complete data structures). */
  virtual void _complete_sweep();

  /*! Initialize an event associated with a point.
   * \param p The given point.
   * \param type The event type.
   */
  void _init_point(const Point_2& pt, Attribute type);

  /*! Initialize the events associated with an x-monotone curve.
   * \param curve The given x-monotone curve.
   * \param index Its unique index.
   */
  void _init_curve(const X_monotone_curve_2& curve, unsigned int index);

  /*! Initialize an event associated with an x-monotone curve end.
   * \param cv The given x-monotone curve.
   * \param ind Its end (ARR_MIN_END or ARR_MAX_END).
   * \param sc The subcurve corresponding to cv.
   */
  void _init_curve_end(const X_monotone_curve_2& cv, Arr_curve_end ind,
                       Subcurve* sc);

  // Variant keeping track of indexed events
  void _init_curve_end(const X_monotone_curve_2& cv, Arr_curve_end ind,
                       Subcurve* sc,
                       std::vector<Event_queue_iterator>& events, std::size_t index);

  /*! Handle the subcurves that are to the left of the event point (i.e.,
   * subcurves that we are done with).
   */
  virtual void _handle_left_curves();

  /*! Handle an event that does not have any incident left curves.
   * Such an event is usually the left endpoint of its incident right
   * subcurves, and we locate their position in the status line.
   */
  void _handle_event_without_left_curves();

  /*! Sort the left subcurves of an event point according to their order in
   * their status line (no geometric comparisons are needed).
   */
  void _sort_left_curves();

  /*! Handle the subcurves to the right of the current event point. */
  virtual void _handle_right_curves();

  /*! Add a subcurve to the right of an event point.
   * \param event The event point.
   * \param curve The subcurve to add.
   * \return (true) if an overlap occurred; (false) otherwise.
   */
  virtual bool _add_curve_to_right(Event* event, Subcurve* curve);

  /*! Add a curve as a right curve or left curve when the event is created
   * or updated.
   */
  virtual void _add_curve(Event* e, Subcurve* sc, Attribute type);

  /*! Remove a curve from the status line. */
  void _remove_curve_from_status_line(Subcurve *leftCurve);

  /*! Allocate an event object associated with a given point.
   * \param pt The point.
   * \param type The event type.
   * \param ps_x The location of the point in x.
   * \param ps_y The location of the point in y.
   * \pre Neither one of the boundary conditions is +/-oo.
   * \return The created event.
   */
  Event* _allocate_event(const Point_2& pt, Attribute type,
                         Arr_parameter_space ps_x, Arr_parameter_space ps_y);

  /*! Allocate an event at open boundary,
   * which is not associated with a valid point.
   * \param type The event type.
   * \param ps_x The location of the point in x.
   * \param ps_y The location of the point in y.
   * \param At least one of the boundary conditions is +/-oo.
   * \return The created event.
   */
  Event* _allocate_event_at_open_boundary(Attribute type,
                                          Arr_parameter_space ps_x,
                                          Arr_parameter_space ps_y);

  /*! Push a finite event point into the event queue.
   * \param pt The point associated with the event.
   * \param type The event type.
   * \param ps_x The location of the point in x.
   * \param ps_y The location of the point in y.
   * \param sc A subcurve that the new event represents on of its endpoints.
   * \return A pair that comprises a pointer to the event, and a flag
   *         indicating whether this is a new event (if false, the event
   *         was in the queue and we just updated it).
   */
  std::pair<Event*, bool> _push_event(const Point_2& pt, Attribute type,
                                      Arr_parameter_space ps_x,
                                      Arr_parameter_space ps_y,
                                      Subcurve* sc = nullptr);

  // Variant keeping track of indexed events
  std::pair<Event*, bool> _push_event(const Point_2& pt, Attribute type,
                                      Arr_parameter_space ps_x,
                                      Arr_parameter_space ps_y,
                                      Subcurve* sc,
                                      std::vector<Event_queue_iterator>& events,
                                      std::size_t index);

  /*! Push an event point associated with a curve end into the event queue.
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
  std::pair<Event*, bool> _push_event(const X_monotone_curve_2& cv,
                                      Arr_curve_end ind,
                                      Attribute type,
                                      Arr_parameter_space ps_x,
                                      Arr_parameter_space ps_y,
                                      Subcurve* sc = nullptr);

  // Variant keeping track of indexed events
  std::pair<Event*, bool> _push_event(const X_monotone_curve_2& cv,
                                      Arr_curve_end ind,
                                      Attribute type,
                                      Arr_parameter_space ps_x,
                                      Arr_parameter_space ps_y,
                                      Subcurve* sc,
                                      const Point_2& pt,
                                      std::vector<Event_queue_iterator>& events,
                                      std::size_t index);

  void _update_event_at_open_boundary(Event* e,
                                      const X_monotone_curve_2& cv,
                                      Arr_curve_end ind,
                                      bool is_new)
  {
    _update_event_at_open_boundary(e, cv, ind, is_new,
                                   Are_all_sides_oblivious_category());
  }

  void _update_event_at_open_boundary(Event* e,
                                      const X_monotone_curve_2& cv,
                                      Arr_curve_end ind,
                                      bool is_new,
                                      Arr_not_all_sides_oblivious_tag)
  { m_visitor->update_event(e, cv, ind, is_new); }

  void _update_event_at_open_boundary(Event* /* e */,
                                      const X_monotone_curve_2& /* cv */,
                                      Arr_curve_end /* ind */,
                                      bool /* is_new */,
                                      Arr_all_sides_oblivious_tag)
  { CGAL_error(); }

#ifdef CGAL_SS_VERBOSE

  uint8_t m_indent_size;
  bool m_need_indent;

  void print_text(const char* text, bool do_eol = false);
  void print_eol();
  void increase_indent();
  void decrease_indent();
  void print_start(const char* name, bool do_eol = true);
  void print_end(const char* name, bool do_eol = true);
  void print_curve(const Subcurve* sc);
  void print_event_info(const Event* e);

  void PrintEventQueue();
  void PrintSubCurves();
  void PrintStatusLine();
  void PrintOpenBoundaryType(Arr_parameter_space x, Arr_parameter_space y);
  void PrintEvent(const Event* e);
#endif

};

} // namespace Surface_sweep_2
} // namespace CGAL

// DEBUG UTILITIES
#ifdef CGAL_SS_VERBOSE
#include <CGAL/Surface_sweep_2/Surface_sweep_2_debug.h>
#endif

#include <CGAL/Surface_sweep_2/No_intersection_surface_sweep_2_impl.h>

#endif
