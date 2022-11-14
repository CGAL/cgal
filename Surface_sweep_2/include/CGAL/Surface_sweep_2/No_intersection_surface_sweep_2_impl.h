// Copyright (c) 2006,2007,2008,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Baruch Zukerman  <baruchzu@post.tau.ac.il>
//                 Efi Fogel        <efif@post.tau.ac.il>
//                 Eric Berberich   <eric.berberich@cgal.org>
//                 (based on old version by Tali Zvi)

#ifndef CGAL_NO_INTERSECTION_SURFACE_SWEEP_2_IMPL_H
#define CGAL_NO_INTERSECTION_SURFACE_SWEEP_2_IMPL_H

#include <CGAL/license/Surface_sweep_2.h>

/*! \file
 *
 * Member-function definitions for the No_intersection_surface_sweep_2 class.
 */

namespace CGAL {
namespace Surface_sweep_2 {

//-----------------------------------------------------------------------------
// Constructor.
//
template <typename Vis>
No_intersection_surface_sweep_2<Vis>::
No_intersection_surface_sweep_2(Visitor* visitor) :
  m_traits(new Traits_adaptor_2()),
  m_traitsOwner(true),
  m_currentEvent(nullptr),
  m_statusLineCurveLess(m_traits, &m_currentEvent),
  m_queueEventLess(m_traits),
  m_queue(new Event_queue(m_queueEventLess)),
  m_statusLine(m_statusLineCurveLess),
  m_status_line_insert_hint(m_statusLine.begin()),
  m_num_of_subCurves(0),
  m_visitor(visitor)
#ifdef CGAL_SS_VERBOSE
  , m_indent_size(0)
  , m_need_indent(false)
#endif
{ m_visitor->attach(this); }

//-----------------------------------------------------------------------------
// Constructor with a given traits-class.
//
template <typename Vis>
No_intersection_surface_sweep_2<Vis>::
No_intersection_surface_sweep_2(const Gt2* traits, Visitor* visitor) :
  m_traits(static_cast<const Traits_adaptor_2*>(traits)),
  m_traitsOwner(false),
  m_currentEvent(nullptr),
  m_statusLineCurveLess(m_traits, &m_currentEvent),
  m_queueEventLess(m_traits),
  m_queue(new Event_queue(m_queueEventLess)),
  m_statusLine(m_statusLineCurveLess),
  m_status_line_insert_hint(m_statusLine.begin()),
  m_num_of_subCurves(0),
  m_visitor(visitor)
#ifdef CGAL_SS_VERBOSE
  , m_indent_size(0)
  , m_need_indent(false)
#endif
{ m_visitor->attach(this); }

//-----------------------------------------------------------------------------
// Destrcutor.
//
template <typename Vis>
No_intersection_surface_sweep_2<Vis>::~No_intersection_surface_sweep_2()
{
  // Free the traits-class object, if we own it.
  if (m_traitsOwner) delete m_traits;

  // Free the event queue.
  delete m_queue;
}

//-----------------------------------------------------------------------------
// Stop the sweep-line process.
//
template <typename Vis>
void No_intersection_surface_sweep_2<Vis>::stop_sweep()
{
  // Clear the event queue, deallocating all events but the current one
  // (the first event in the queue).
  Event_queue_iterator qiter= this->m_queue->begin();

  for (++qiter; qiter != this->m_queue->end(); ++qiter)
    this->deallocate_event(*qiter);

  // Clear the status line.
  this-> m_statusLine.clear();
  m_status_line_insert_hint = this->m_statusLine.begin();

  // Empty the event queue, and leave only the first event there.
  CGAL_assertion(!m_queue->empty());
  Event_queue_iterator second = m_queue->begin();
  ++second;
  while (second != m_queue->end()) {
    Event_queue_iterator next = second;
    ++next;
    m_queue->erase(second);
    second = next;
  }
}

//-----------------------------------------------------------------------------
// Deallocate event object..
//
template <typename Vis>
void No_intersection_surface_sweep_2<Vis>::deallocate_event(Event* event)
{
  // Remove the event from the set of allocated events.
  m_allocated_events.erase(m_allocated_events.iterator_to(*event));
}

//-----------------------------------------------------------------------------
// Perform the main sweep-line loop.
//
template <typename Vis>
void No_intersection_surface_sweep_2<Vis>::_sweep()
{
  CGAL_SS_PRINT_TEXT("Ordered sequence of ");
  CGAL_SS_PRINT(m_queue->size());
  CGAL_SS_PRINT_TEXT(" initial events:");
  CGAL_SS_PRINT_EOL();
  CGAL_SS_DEBUG(Event_queue_iterator eventIter1 = m_queue->begin();
                while (eventIter1 != m_queue->end()) {
                  print_text("* ");
                  PrintEvent(*eventIter1);
                  print_eol();
                  eventIter1++;
                })

  // Looping over the events in the queue.
  Event_queue_iterator eventIter = m_queue->begin();

  while (eventIter != m_queue->end()) {
    // Get the next event from the queue.
    m_currentEvent = *eventIter;

    CGAL_SS_PRINT_EOL();
    CGAL_SS_PRINT_TEXT("------------- ");
    CGAL_SS_DEBUG(PrintEvent(m_currentEvent));
    CGAL_SS_PRINT_TEXT( " --------------");
    CGAL_SS_PRINT_EOL();
    CGAL_SS_PRINT_EOL();
    CGAL_SS_PRINT_STATUS_LINE();
    CGAL_SS_PRINT_EVENT_INFO(m_currentEvent);

    // Handle the subcurves that are to the left of the event point (i.e.,
    // subcurves that we are done with).
    _handle_left_curves();

    // Handle the subcurves to the right of the event point, reorder them
    // and test for intersections between them and their immediate neighbors
    // on the status line.
    _handle_right_curves();

    // Inform the visitor about the event. The visitor also determines whether
    // it is possible to deallocate the event now, or will it be deallocated
    // later (at the visitor's responsibility).
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
}

//-----------------------------------------------------------------------------
// Initialize the data structures for the sweep-line algorithm.
//
template <typename Vis>
void No_intersection_surface_sweep_2<Vis>::_init_structures()
{
  CGAL_assertion(m_queue->empty());
  CGAL_assertion((m_statusLine.size() == 0));

  // Allocate all of the Subcurve objects as one block. Don't allocate
  // anything when there are no subcurves.
  if (m_num_of_subCurves > 0)
    m_subCurves = m_subCurveAlloc.allocate(m_num_of_subCurves);
}

//-----------------------------------------------------------------------------
// Complete the sweep (complete the data structures).
//
template <typename Vis>
void No_intersection_surface_sweep_2<Vis>::_complete_sweep()
{
  CGAL_assertion(m_queue->empty());
  CGAL_assertion((m_statusLine.size() == 0));

  // Free all subcurve objects.
  for (unsigned int i = 0; i < m_num_of_subCurves; ++i){
    std::allocator_traits<Subcurve_alloc>::destroy(m_subCurveAlloc, m_subCurves + i);
  }

  if (m_num_of_subCurves > 0)
    m_subCurveAlloc.deallocate(m_subCurves, m_num_of_subCurves);
}

//-----------------------------------------------------------------------------
// Initialize an event associated with a point.
//
template <typename Vis>
void No_intersection_surface_sweep_2<Vis>::_init_point(const Point_2& pt,
                                                       Attribute type)
{
  // Create the event, or obtain an existing event in the queue.
  Arr_parameter_space ps_x = m_traits->parameter_space_in_x_2_object()(pt);
  Arr_parameter_space ps_y = m_traits->parameter_space_in_y_2_object()(pt);
#if 0
  CGAL::set_pretty_mode(std::cout);
  std::cout << "init pt ps_x: " << ps_x << std::endl;
  std::cout << "init pt ps_y: " << ps_y << std::endl;
#endif

  const std::pair<Event*, bool>& pair_res =
    _push_event(pt, type, ps_x, ps_y);

  bool is_new = pair_res.second;
  m_visitor->update_event(pair_res.first, pt, is_new);
}

//-----------------------------------------------------------------------------
// Initialize the events associated with an x-monotone curve.
//
template <typename Vis>
void No_intersection_surface_sweep_2<Vis>::
_init_curve(const X_monotone_curve_2& curve, unsigned int index)
{
  // Construct and initialize a subcurve object.
  Subcurve* sc = m_subCurves + index;
  std::allocator_traits<Subcurve_alloc>::construct(m_subCurveAlloc, sc, m_masterSubcurve );
  sc->set_hint(this->m_statusLine.end());
  sc->init(curve);

  // Create two events associated with the curve ends, respectively.
  _init_curve_end(curve, ARR_MAX_END, sc, All_sides_oblivious_category());
  _init_curve_end(curve, ARR_MIN_END, sc, All_sides_oblivious_category());
}

//-----------------------------------------------------------------------------
// Initialize an event associated with an x-monotone curve end.
// This is the implementation for the case where all 4 boundary sides are
// oblivious.
//
template <typename Vis>
void No_intersection_surface_sweep_2<Vis>::
_init_curve_end(const X_monotone_curve_2& cv, Arr_curve_end ind, Subcurve* sc,
                Arr_all_sides_oblivious_tag)
{
  const Attribute end_attr =
    (ind == ARR_MIN_END) ? Event::LEFT_END : Event::RIGHT_END;

  const Point_2& pt = (ind == ARR_MIN_END) ?
    m_traits->construct_min_vertex_2_object()(cv) :
    m_traits->construct_max_vertex_2_object()(cv);

  // Create the corresponding event and push it into the event queue.
  std::pair<Event*, bool> pair_res =
    _push_event(pt, end_attr, ARR_INTERIOR, ARR_INTERIOR, sc);

  // Inform the visitor in case we updated an existing event.
  m_visitor->update_event(pair_res.first, pt, cv, ind, pair_res.second);
}

//-----------------------------------------------------------------------------
// Initialize an event associated with an x-monotone curve end.
// This is the implementation for the case where there is at least one boundary
// side that is not oblivious.
//
//
template <typename Vis>
void No_intersection_surface_sweep_2<Vis>::
_init_curve_end(const X_monotone_curve_2& cv, Arr_curve_end ind, Subcurve* sc,
                Arr_not_all_sides_oblivious_tag)
{
  const Attribute end_attr =
    (ind == ARR_MIN_END) ? Event::LEFT_END : Event::RIGHT_END;

  // Get the parameter space of the curve end.
  Arr_parameter_space ps_x = m_traits->parameter_space_in_x_2_object()(cv, ind);
  Arr_parameter_space ps_y = m_traits->parameter_space_in_y_2_object()(cv, ind);
  // Create the corresponding event and push it into the event queue.

  if (m_traits->is_closed_2_object()(cv, ind)) {
    // The curve end is closed and thus associated with a valid endpoint.
    const Point_2& pt = (ind == ARR_MIN_END) ?
      m_traits->construct_min_vertex_2_object()(cv) :
      m_traits->construct_max_vertex_2_object()(cv);

    // Create the corresponding event and push it into the event queue.
    std::pair<Event*, bool> pair_res =
      ((ps_x == ARR_INTERIOR) && (ps_y == ARR_INTERIOR)) ?
      _push_event(pt, end_attr, ps_x, ps_y, sc) :
      _push_event(cv, ind, end_attr, ps_x, ps_y, sc);

    // Inform the visitor in case we updated an existing event.
    Event* e = pair_res.first;
    CGAL_assertion(e->is_closed());
    m_visitor->update_event(e, pt, cv, ind, pair_res.second);
    return;
  }

  // The curve end is open.
  // Create the corresponding event and push it into the event queue.
  std::pair<Event*, bool> pair_res =
    _push_event(cv, ind, end_attr, ps_x, ps_y, sc);

  // Inform the visitor in case we updated an existing event.
  Event* e = pair_res.first;
  CGAL_assertion(! e->is_closed());
  m_visitor->update_event(e, cv, ind, pair_res.second);
}

template <typename Vis>
void No_intersection_surface_sweep_2<Vis>::
_init_curve_end(const X_monotone_curve_2& cv, Arr_curve_end ind, Subcurve* sc,
                std::vector<Event_queue_iterator>& events, std::size_t index)
{
  // Get the boundary conditions of the curve end.
  const Attribute  end_attr =
    (ind == ARR_MIN_END) ? Event::LEFT_END : Event::RIGHT_END;

  Arr_parameter_space ps_x = m_traits->parameter_space_in_x_2_object()(cv, ind);
  Arr_parameter_space ps_y = m_traits->parameter_space_in_y_2_object()(cv, ind);

  // Create the corresponding event and push it into the event queue.
  std::pair<Event*, bool> pair_res;

  if (m_traits->is_closed_2_object()(cv, ind)) {
    // The curve end is closed and thus associated with a valid endpoint.
    const Point_2& pt = (ind == ARR_MIN_END) ?
      m_traits->construct_min_vertex_2_object()(cv) :
      m_traits->construct_max_vertex_2_object()(cv);

    pair_res = ((ps_x == ARR_INTERIOR) && (ps_y == ARR_INTERIOR)) ?
      _push_event(pt, end_attr, ps_x, ps_y, sc, events, index) :
      _push_event(cv, ind, end_attr, ps_x, ps_y, sc, pt, events, index);

    // Inform the visitor in case we updated an existing event.
    Event* e = pair_res.first;
    CGAL_assertion(e->is_closed());
    m_visitor->update_event(e, pt, cv, ind, pair_res.second);
  }
  else {
    // The curve end is open, insert it into the event queue.
    pair_res = _push_event(cv, ind, end_attr, ps_x, ps_y, sc);

    // Inform the visitor in case we updated an existing event.
    Event* e = pair_res.first;
    CGAL_assertion(! e->is_closed());
    _update_event_at_open_boundary(e, cv, ind, pair_res.second);
  }
}

//-----------------------------------------------------------------------------
// Handle the subcurves to the left of the current event point.
//
template <typename Vis>
void No_intersection_surface_sweep_2<Vis>::_handle_left_curves()
{
  CGAL_SS_PRINT_START_EOL("handling left curves at (");
  CGAL_SS_DEBUG(PrintEvent(m_currentEvent));
  CGAL_SS_PRINT_EOL();

  m_is_event_on_above = false;

  if (! m_currentEvent->has_left_curves()) {
    // The current event does not have any incident left subcurves, so we
    // should find a place for it in the status line (the function we call
    // update the m_status_line_insert_hint and m_is_event_on_above members).
    // We also notify the visitor on the new event we are about to handle.
    _handle_event_without_left_curves(Sides_category());
    m_visitor->before_handle_event(m_currentEvent);

    // Nothing else to do (no left curves).
    CGAL_SS_PRINT_END_EOL("handling left curves");
    return;
  }

  CGAL_SS_PRINT_TEXT("left curves before sorting:");
  CGAL_SS_PRINT_EOL();
  CGAL_SS_DEBUG(if (m_currentEvent->left_curves_begin() !=
                    m_currentEvent->left_curves_end())
                { print_event_info(m_currentEvent); });

  // Use the status-line to sort all left subcurves incident to the current
  // event (no geometric comparisons are neede at all).
  _sort_left_curves();

  // Now the event is updated, with its left subcurved properly sorted, and
  // we can inform the visitor that we are about to handle this event.
  m_visitor->before_handle_event(m_currentEvent);

  CGAL_SS_PRINT_TEXT("left curves after sorting:");
  CGAL_SS_PRINT_EOL();
  CGAL_SS_DEBUG(if (m_currentEvent->left_curves_begin() !=
                    m_currentEvent->left_curves_end() )
                { print_event_info(m_currentEvent); });

  // Remove all left subcurves from the status line, and inform the visitor
  // that we are done handling these subcurves.
  Event_subcurve_iterator left_iter = m_currentEvent->left_curves_begin();

  while (left_iter != m_currentEvent->left_curves_end()) {
    Subcurve* left_sc = *left_iter;

    m_visitor->add_subcurve(left_sc->last_curve(), left_sc);
    ++left_iter;

    _remove_curve_from_status_line(left_sc);
  }
  CGAL_SS_PRINT_END_EOL("handling left curves");
}

//-----------------------------------------------------------------------------
// Handle an event that does not have any incident left curves.
//
template <typename Vis>
void No_intersection_surface_sweep_2<Vis>::
_handle_event_without_left_curves(Arr_all_sides_oblivious_tag)
{
  const std::pair<Status_line_iterator, bool>& pair_res =
    m_statusLine.find_lower(m_currentEvent->point(), m_statusLineCurveLess);
  m_status_line_insert_hint = pair_res.first;
  m_is_event_on_above = pair_res.second;
}

//-----------------------------------------------------------------------------
// Handle an event that does not have any incident left curves.
//
template <typename Vis>
void No_intersection_surface_sweep_2<Vis>::
_handle_event_without_left_curves(Arr_all_sides_not_finite_tag)
{
  // Check if the event is a boundary event or not.
  const Arr_parameter_space ps_x = m_currentEvent->parameter_space_in_x();
  const Arr_parameter_space ps_y = m_currentEvent->parameter_space_in_y();

  if ((ps_x == ARR_INTERIOR) && (ps_y == ARR_INTERIOR)) {
    _handle_event_without_left_curves(Arr_all_sides_oblivious_tag());
    return;
  }

  // Process a boundary event.
  if ((ps_x == ARR_LEFT_BOUNDARY) || (ps_y == ARR_TOP_BOUNDARY)) {
    // If the event is on the left boundary, then we are still sweeping the
    // left boundary. While we sweep the left boundary we maintain the order
    // of the events in the queue; thus, the new event should be placed above
    // all other subcurves in the status line.
    // If the event is on the top, it should also be inserted below all
    // other subcurves
    m_status_line_insert_hint = m_statusLine.end();
    return;
  }

  // Note that an event on the right boundary can only represent a right end
  // of a curve.
  // If the event is on the bottom boundary, it should be inserted below all
  // other subcurves;
  CGAL_assertion(ps_x != ARR_RIGHT_BOUNDARY);
  CGAL_assertion(ps_y == ARR_BOTTOM_BOUNDARY);
  m_status_line_insert_hint = m_statusLine.begin();
}

//-----------------------------------------------------------------------------
// Handle an event that does not have any incident left curves.
//
template <typename Vis>
void No_intersection_surface_sweep_2<Vis>::
_handle_event_without_left_curves(Arr_not_all_sides_not_finite_tag)
{
  // Check if the event is a boundary event or not.
  const Arr_parameter_space ps_x = m_currentEvent->parameter_space_in_x();
  const Arr_parameter_space ps_y = m_currentEvent->parameter_space_in_y();

  if ((ps_x == ARR_INTERIOR) && (ps_y == ARR_INTERIOR)) {
    _handle_event_without_left_curves(Arr_all_sides_oblivious_tag());
    return;
  }

  // Process a boundary event.
  if (ps_x == ARR_LEFT_BOUNDARY) {
    // If the event is on the left boundary, then we are still sweeping the
    // left boundary. While we sweep the left boundary we maintain the order
    // of the events in the queue; thus, the new event should be placed above
    // all other subcurves in the status line, unless the status line
    // contains a curve that entirely lies on the left boundary.
    if (m_currentEvent->is_closed()) {
      _handle_event_without_left_curves(Arr_all_sides_oblivious_tag());
      return;
    }
    m_status_line_insert_hint = m_statusLine.end();
    return;
  }

  // Note that an event with a positive boundary condition at x can only
  // represent a right end of a curve.
  CGAL_assertion(ps_x != ARR_RIGHT_BOUNDARY);

  // If the event is on the bottom boundary, it should be inserted below all
  // other subcurves.
  if (ps_y == ARR_BOTTOM_BOUNDARY) {
    m_status_line_insert_hint = m_statusLine.begin();
    return;
  }

  // If the event is on the top boundary, it should be inserted above all
  // other subcurves.
  CGAL_assertion(ps_y == ARR_TOP_BOUNDARY);
  m_status_line_insert_hint = m_statusLine.end();
}

//-----------------------------------------------------------------------------
// Sort the left subcurves of an event point according to their order in
// their status line (no geometric comprasions are needed).
//
template <typename Vis>
void No_intersection_surface_sweep_2<Vis>::_sort_left_curves()
{
  CGAL_SS_PRINT_START_EOL("sorting left curves");
  CGAL_assertion(m_currentEvent->has_left_curves());

  // Get the first curve associated with the event and its position on the
  // status line. We proceed from this position up the status line until
  // we encounter a subcurve that is not associated with the current event.
  Subcurve* curve = *(m_currentEvent->left_curves_begin());
  Status_line_iterator sl_iter = curve->hint();

  CGAL_assertion(sl_iter != m_statusLine.end());
  CGAL_assertion(*sl_iter == curve);
  // Look for the first curve in the vertical ordering that is also in the
  // left curve of the event
  Status_line_iterator end = sl_iter;
  for (++end; end != m_statusLine.end(); ++end) {
    if (std::find(m_currentEvent->left_curves_begin(),
                  m_currentEvent->left_curves_end(), *end) ==
        m_currentEvent->left_curves_end())
      break;
  }

  if (sl_iter == m_statusLine.begin()) {
    // In case the lowest subcurve in the status line is associated with the
    // current event, we have the range of (sorted) subcurves ready. We
    // associate this range with the event, so the curves are now sorted
    // according to their vertical positions immediately to the left of the
    // event.
    m_currentEvent->replace_left_curves(sl_iter, end);
    CGAL_SS_PRINT_END_EOL("sorting left curves");
    return;
  }

  // Go down the status line until we encounter a subcurve that is not
  // associated with the current event.
  for (--sl_iter; sl_iter != m_statusLine.begin(); --sl_iter) {
    if (std::find(m_currentEvent->left_curves_begin(),
                  m_currentEvent->left_curves_end(), *sl_iter) ==
        m_currentEvent->left_curves_end())
    {
      // Associate the sorted range of subcurves with the event.
      m_currentEvent->replace_left_curves(++sl_iter, end);
      CGAL_SS_PRINT_END_EOL("sorting left curves");
      return;
    }
  }

  // Check if the subcurve at the current iterator position should be
  // associated with the current event, and select the (sorted) range of
  // subcurves accordingly.
  if (std::find(m_currentEvent->left_curves_begin(),
                m_currentEvent->left_curves_end(), *sl_iter) ==
      m_currentEvent->left_curves_end())
    m_currentEvent->replace_left_curves(++sl_iter, end);
  else
    m_currentEvent->replace_left_curves(sl_iter, end);

  CGAL_SS_PRINT_END_EOL("sorting left curves");
}

//-----------------------------------------------------------------------------
// Handle the subcurves to the right of the current event point.
//
template <typename Vis>
void No_intersection_surface_sweep_2<Vis>::_handle_right_curves()
{
  CGAL_SS_PRINT_START("Handling right curves at (");
  CGAL_SS_DEBUG(PrintEvent(m_currentEvent));
  CGAL_SS_PRINT_EOL();

  // We have nothing to do if the current event does not have any incident
  // right subcurves.
  if (! m_currentEvent->has_right_curves()) return;

  // Loop over the curves to the right of the current event and handle them:
  // Since the curves are non intersecting, the event can represents the
  // left end of the right curves and we have no prior information from the
  // order of the left subcurves. Thus, we just insert the curves to the
  // status line.
  Event_subcurve_iterator curr = m_currentEvent->right_curves_begin();
  Event_subcurve_iterator right_end = m_currentEvent->right_curves_end();

  while (curr != right_end) {
    Subcurve* sc = *curr;
    CGAL_SS_PRINT_INSERT(sc);

    // Insert the curve to the left-curves of the right event.
    // sc->right_event()->add_curve_to_left(sc);

    // Insert the curve into the status line.
    Status_line_iterator sl_iter =
      m_statusLine.insert_before(m_status_line_insert_hint, sc);
    sc->set_hint(sl_iter);

    CGAL_SS_DEBUG(PrintStatusLine(););
    ++curr;
  }

  CGAL_SS_PRINT_STATUS_LINE();
  CGAL_SS_PRINT_END_EOL("handling right curves done");
}

//-----------------------------------------------------------------------------
// Add a subcurve to the right of an event point.
//
template <typename Vis>
bool No_intersection_surface_sweep_2<Vis>::_add_curve_to_right(Event* event,
                                                               Subcurve* curve)
{
#if defined(CGAL_NO_ASSERTIONS)
  (void) event->add_curve_to_right(curve, m_traits);
#else
  std::pair<bool, Event_subcurve_iterator> pair_res =
     event->add_curve_to_right(curve, m_traits);
  CGAL_assertion(!pair_res.first);
#endif

  return false;
}

//-----------------------------------------------------------------------------
// Remove a curve from the status line.
//
template <typename Vis>
void No_intersection_surface_sweep_2<Vis>::
_remove_curve_from_status_line(Subcurve* sc)
{
  CGAL_SS_PRINT_START("removing a curve from the status line, ");
  CGAL_SS_PRINT_CURVE(sc);
  CGAL_SS_PRINT_EOL();
  CGAL_SS_PRINT_STATUS_LINE();

  // Get the position of the subcurve on the status line.
  Status_line_iterator sl_iter = sc->hint();
  CGAL_assertion(sl_iter != m_statusLine.end());

  // The position of the next event can be right after the deleted subcurve.
  m_status_line_insert_hint = sl_iter;
  ++m_status_line_insert_hint;
  sc->set_hint(m_statusLine.end());

  // Erase the subcurve from the status line.
  CGAL_SS_PRINT_ERASE(*sl_iter);
  m_statusLine.erase(sl_iter);
  CGAL_SS_PRINT_END_EOL("removing a curve from the status line");
}

//-----------------------------------------------------------------------------
// Allocate an event object associated with a valid point.
//
template <typename Vis>
typename No_intersection_surface_sweep_2<Vis>::Event*
No_intersection_surface_sweep_2<Vis>::_allocate_event(const Point_2& pt,
                                                      Attribute type,
                                                      Arr_parameter_space ps_x,
                                                      Arr_parameter_space ps_y)
{
  // Allocate the event.
  Event* e = &*m_allocated_events.emplace();
  e->init(pt, type, ps_x, ps_y);
  return e;
}

//-----------------------------------------------------------------------------
// Allocate an event at open boundary,
// which is not associated with a valid point.
//
template <typename Vis>
typename No_intersection_surface_sweep_2<Vis>::Event*
No_intersection_surface_sweep_2<Vis>::
_allocate_event_at_open_boundary(Attribute type,
                                 Arr_parameter_space ps_x,
                                 Arr_parameter_space ps_y)
{
  Event* e = &*m_allocated_events.emplace();
  e->init_at_open_boundary(type, ps_x, ps_y);
  return e;
}

//-----------------------------------------------------------------------------
// Push a closed event point into the event queue.
//
template <typename Vis>
std::pair<typename No_intersection_surface_sweep_2<Vis>::Event*, bool>
No_intersection_surface_sweep_2<Vis>::_push_event(const Point_2& pt,
                                                  Attribute type,
                                                  Arr_parameter_space ps_x,
                                                  Arr_parameter_space ps_y,
                                                  Subcurve* sc)
{
  // Look for the point in the event queue.
  Event* e;
  m_queueEventLess.set_parameter_space_in_x(ps_x);
  m_queueEventLess.set_parameter_space_in_y(ps_y);

  const std::pair<Event_queue_iterator, bool>& pair_res =
    m_queue->find_lower(pt, m_queueEventLess);
  const bool exist = pair_res.second;
  if (!exist) {
    // The point is not found in the event queue - create a new event and
    // insert it into the queue.
    e = _allocate_event(pt, type, ps_x, ps_y);
  }
  else {
    // The event associated with the given point already exists in the queue,
    // so we just have to update it.
    e = *(pair_res.first);
    CGAL_assertion(e->is_closed());

    e->set_attribute(type);
  }
  CGAL_assertion(e->parameter_space_in_x() == ps_x);
  CGAL_assertion(e->parameter_space_in_y() == ps_y);

  // If we are given a subcurve that the event represents one of its
  // endpoints, update the event and the subcurve records accordingly.
  // Note that this must be done before we actually insert the new event
  // into the event queue.
  _add_curve(e, sc, type);

  // Insert the new event into the queue using the hint we got when we
  // looked for it.
  if (! exist) m_queue->insert_before(pair_res.first, e);

#ifdef CGAL_SS_VERBOSE
  if (! exist) {
    CGAL_SS_PRINT_NEW_EVENT(pt, e);
  }
  else {
    CGAL_SS_PRINT_UPDATE_EVENT(pt, e);
  }
#endif

  // Return the resulting event and a flag indicating whether we have created
  // a new event.
  return (std::make_pair(e, !exist));
}

template <typename Vis>
std::pair<typename No_intersection_surface_sweep_2<Vis>::Event*, bool>
No_intersection_surface_sweep_2<Vis>::_push_event(const Point_2& pt,
                                                  Attribute type,
                                                  Arr_parameter_space ps_x,
                                                  Arr_parameter_space ps_y,
                                                  Subcurve* sc,
                                                  std::vector<Event_queue_iterator>& events,
                                                  std::size_t index)
{
  Event* e;

  std::pair<Event_queue_iterator, bool>
    pair_res = std::make_pair (events[index], true);

  // If event does not exist
  if (events[index] == Event_queue_iterator())
  {
    // Still look for the curve end in the event queue in case two
    // point are the the same in the vertex range
    m_queueEventLess.set_parameter_space_in_x(ps_x);
    m_queueEventLess.set_parameter_space_in_y(ps_y);
    pair_res = m_queue->find_lower(pt, m_queueEventLess);
  }

  bool exist = pair_res.second;
  if (! exist) {
    // The point is not found in the event queue - create a new event and
    // insert it into the queue.
    e = _allocate_event(pt, type, ps_x, ps_y);
  }
  else {
    events[index] = pair_res.first;
    // The event associated with the given point already exists in the queue,
    // so we just have to update it.
    e = *(pair_res.first);
    CGAL_assertion(e->is_closed());

    e->set_attribute(type);
  }

  // If we are given a subcurve that the event represents one of its
  // endpoints, update the event and the subcurve records accordingly.
  // Note that this must be done before we actually insert the new event
  // into the event queue.
  _add_curve(e, sc, type);

  // Insert the new event into the queue using the hint we got when we
  // looked for it.
  if (! exist)
    events[index] = m_queue->insert_before(pair_res.first, e);

#ifdef CGAL_SS_VERBOSE
  if (! exist) {
    CGAL_SS_PRINT_NEW_EVENT(pt, e);
  }
  else {
    CGAL_SS_PRINT_UPDATE_EVENT(pt, e);
  }
#endif

  // Return the resulting event and a flag indicating whether we have created
  // a new event.
  return (std::make_pair(e, !exist));
}

//-----------------------------------------------------------------------------
// Push an event point associated with a curve end into the event queue.
//
template <typename Vis>
std::pair<typename No_intersection_surface_sweep_2<Vis>::Event*, bool>
No_intersection_surface_sweep_2<Vis>::_push_event(const X_monotone_curve_2& cv,
                                                  Arr_curve_end ind,
                                                  Attribute type,
                                                  Arr_parameter_space ps_x,
                                                  Arr_parameter_space ps_y,
                                                  Subcurve* sc)
{
  // Look for the curve end in the event queue.
  Event* e;

  m_queueEventLess.set_parameter_space_in_x(ps_x);
  m_queueEventLess.set_parameter_space_in_y(ps_y);
  m_queueEventLess.set_index(ind);

  const std::pair<Event_queue_iterator, bool>& pair_res =
    m_queue->find_lower(cv, m_queueEventLess);
  const bool exist = pair_res.second;

  if (! exist) {
    // The curve end is not found in the event queue - create a new event and
    // insert it into the queue.
    if (m_traits->is_closed_2_object()(cv, ind)) {
      // The curve end is closed and so it is associated with a valid
      // point.
      const Point_2& pt = (ind == ARR_MIN_END) ?
        m_traits->construct_min_vertex_2_object()(cv) :
        m_traits->construct_max_vertex_2_object()(cv);

      e = _allocate_event(pt, type, ps_x, ps_y);
    }
    else {
      // The curve end is open, so we create an event at open boundary.
      e = _allocate_event_at_open_boundary(type, ps_x, ps_y);
    }
  }
  else {
    // The event associated with the given curve end already exists in the
    // queue, so we just have to update it.
    e = *(pair_res.first);
#if 0
    std::cout << "ps_x: " << ps_x << std::endl;
    std::cout << "ps_y: " << ps_y << std::endl;
    std::cout << "es_x: " << e->parameter_space_in_x() << std::endl;
    std::cout << "es_y: " << e->parameter_space_in_y() << std::endl;
#endif
    CGAL_assertion(e->parameter_space_in_x() == ps_x);
    CGAL_assertion(e->parameter_space_in_y() == ps_y);

    e->set_attribute(type);
  }

  // If we are given a subcurve that the event represents one of its
  // endpoints, update the event and the subcurve records accordingly.
  // Note that this must be done before we actually insert the new event
  // into the event queue.
  _add_curve(e, sc, type);

  // Insert the new event into the queue using the hint we got when we
  // looked for it.
  if (! exist) m_queue->insert_before(pair_res.first, e);

  return (std::make_pair(e, !exist));
}

template <typename Vis>
std::pair<typename No_intersection_surface_sweep_2<Vis>::Event*, bool>
No_intersection_surface_sweep_2<Vis>::_push_event(const X_monotone_curve_2& cv,
                                                  Arr_curve_end ind,
                                                  Attribute type,
                                                  Arr_parameter_space ps_x,
                                                  Arr_parameter_space ps_y,
                                                  Subcurve* sc,
                                                  const Point_2& pt,
                                                  std::vector<Event_queue_iterator>& events,
                                                  std::size_t index)
{
  Event* e;

  std::pair<Event_queue_iterator, bool>
    pair_res = std::make_pair (events[index], true);

  // If event does not exist
  if (events[index] == Event_queue_iterator())
  {
    // Still look for the curve end in the event queue in case two
    // point are the the same in the vertex range

    m_queueEventLess.set_parameter_space_in_x(ps_x);
    m_queueEventLess.set_parameter_space_in_y(ps_y);
    m_queueEventLess.set_index(ind);

    pair_res =
      m_queue->find_lower(cv, m_queueEventLess);
  }

  bool exist = pair_res.second;

  if (! exist) {
    // The curve end is not found in the event queue - create a new event and
    // insert it into the queue.
    // The curve end is closed and so it is associated with a valid
    // point.
    e = _allocate_event(pt, type, ps_x, ps_y);
  }
  else {
    events[index] = pair_res.first;

    // The event associated with the given curve end already exists in the
    // queue, so we just have to update it.
    e = *(pair_res.first);
    CGAL_assertion((e->parameter_space_in_x() == ps_x) &&
                   (e->parameter_space_in_y() == ps_y));

    e->set_attribute(type);
  }

  // If we are given a subcurve that the event represents one of its
  // endpoints, update the event and the subcurve records accordingly.
  // Note that this must be done before we actually insert the new event
  // into the event queue.
  _add_curve(e, sc, type);

  // Insert the new event into the queue using the hint we got when we
  // looked for it.
  if (! exist)
    events[index] = m_queue->insert_before(pair_res.first, e);

  return (std::make_pair(e, !exist));
}

//-----------------------------------------------------------------------------
// add a curve as a right curve or left curve when the event is created
// or updated.
//
template <typename Vis>
void No_intersection_surface_sweep_2<Vis>::_add_curve(Event* e, Subcurve* sc,
                                                      Attribute type)
{
  if (sc == nullptr) return;

  if (type == Event::LEFT_END) {
    sc->set_left_event(e);
    _add_curve_to_right(e, sc);
    return;
  }

  CGAL_assertion(type == Event::RIGHT_END);
  sc->set_right_event(e);
  // Defer the insertion of the curve to the left-curves of the right event.
  e->add_curve_to_left(sc);
}

} // namespace Surface_sweep_2
} // namespace CGAL

#endif
