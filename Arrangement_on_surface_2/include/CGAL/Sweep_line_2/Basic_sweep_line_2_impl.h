// Copyright (c) 2006,2007,2008,2009,2010,2011 Tel-Aviv University (Israel).
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
// Author(s)     : Baruch Zukerman  <baruchzu@post.tau.ac.il>
//                 Efi Fogel        <efif@post.tau.ac.il>
//                 Eric Berberich   <ericb@post.tau.ac.il>
//                 (based on old version by Tali Zvi)

#ifndef CGAL_BASIC_SWEEP_LINE_2_IMPL_H
#define CGAL_BASIC_SWEEP_LINE_2_IMPL_H

/*! \file
 * Member-function definitions for the Basic_sweep_line_2 class.
 */

namespace CGAL {

//-----------------------------------------------------------------------------
// Constructor.
//
template <typename Tr, typename Vis, typename Subcv, typename Evnt,
          typename Alloc>
Basic_sweep_line_2<Tr, Vis, Subcv, Evnt, Alloc>::
Basic_sweep_line_2(Visitor* visitor) :
  m_traits(new Traits_adaptor_2()),
  m_traitsOwner(true),
  m_statusLineCurveLess(m_traits, &m_currentEvent),
  m_queueEventLess(m_traits),
  m_queue(new Event_queue(m_queueEventLess)),
  m_statusLine(m_statusLineCurveLess),
  m_status_line_insert_hint(m_statusLine.begin()),
  m_num_of_subCurves(0),
  m_visitor(visitor)
{
  m_visitor->attach(this);
}

//-----------------------------------------------------------------------------
// Constructor with a given traits-class.
//
template <typename Tr, typename Vis, typename Subcv, typename Evnt,
          typename Alloc>
Basic_sweep_line_2<Tr, Vis, Subcv, Evnt, Alloc>::
Basic_sweep_line_2(const Traits_2* traits, Visitor* visitor) :
  m_traits(static_cast<const Traits_adaptor_2*>(traits)),
  m_traitsOwner(false),
  m_statusLineCurveLess(m_traits, &m_currentEvent),
  m_queueEventLess(m_traits),
  m_queue(new Event_queue(m_queueEventLess)),
  m_statusLine(m_statusLineCurveLess),
  m_status_line_insert_hint(m_statusLine.begin()),
  m_num_of_subCurves(0),
  m_visitor(visitor)
{
  m_visitor->attach(this);
}

//-----------------------------------------------------------------------------
// Destrcutor.
//
template <typename Tr, typename Vis, typename Subcv, typename Evnt,
          typename Alloc>
Basic_sweep_line_2<Tr, Vis, Subcv, Evnt, Alloc>::~Basic_sweep_line_2()
{
  // Free the traits-class object, if we own it.
  if (m_traitsOwner) delete m_traits;

  // Free the event queue.
  delete m_queue;

  // Free all the event that have not been de-allocated so far.
  Allocated_events_iterator iter;
  Event* p_event;
  for (iter = m_allocated_events.begin();
       iter != m_allocated_events.end(); ++iter)
  {
    p_event = *iter;
    m_eventAlloc.destroy(p_event);
    m_eventAlloc.deallocate(p_event,1);
  }
}

//-----------------------------------------------------------------------------
// Stop the sweep-line process.
//
template <typename Tr, typename Vis, typename Subcv, typename Evnt,
          typename Alloc>
void Basic_sweep_line_2<Tr, Vis, Subcv, Evnt, Alloc>::stop_sweep()
{
  // Clear the event queue, deallocating all events but the current one
  // (the first event in the queue).
  Event_queue_iterator qiter= this->m_queue->begin();

  for (++qiter; qiter != this->m_queue->end(); ++qiter)
    this ->deallocate_event(*qiter);

  // Clear the status line.
  this -> m_statusLine.clear();
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
template <typename Tr, typename Vis, typename Subcv, typename Evnt,
          typename Alloc>
void Basic_sweep_line_2<Tr, Vis, Subcv, Evnt, Alloc>::
deallocate_event(Event* event)
{
  // Remove the event from the set of allocated events.
  m_allocated_events.erase(event);
  
  // Perfrom the actual deallocation.
  m_eventAlloc.destroy(event);
  m_eventAlloc.deallocate(event,1);
}

//-----------------------------------------------------------------------------
// Perform the main sweep-line loop.
//
template <typename Tr, typename Vis, typename Subcv, typename Evnt,
          typename Alloc>
void Basic_sweep_line_2<Tr, Vis, Subcv, Evnt, Alloc>::_sweep()
{
  CGAL_SL_DEBUG(
    {
      CGAL_PRINT("Ordered sequence of " << m_queue->size() 
                 <<  " initial events:\n");
      Event_queue_iterator eventIter1 = m_queue->begin();
      while (eventIter1 != m_queue->end()) {
          
        CGAL_PRINT("* ");
        CGAL_SL_DEBUG(PrintEvent(*eventIter1););
        CGAL_PRINT( "\n");
        eventIter1++;
      }
    }
  )

  // Looping over the events in the queue.
  Event_queue_iterator eventIter = m_queue->begin();
  
  while (eventIter != m_queue->end()) {
    // Get the next event from the queue.
    m_currentEvent = *eventIter;
    
    CGAL_PRINT("------------- ");
    CGAL_SL_DEBUG(PrintEvent(m_currentEvent););
    CGAL_PRINT( " --------------\n");
    CGAL_SL_DEBUG(PrintStatusLine();
                  m_currentEvent->Print(););
      
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
template <typename Tr, typename Vis, typename Subcv, typename Evnt,
          typename Alloc>
void Basic_sweep_line_2<Tr, Vis, Subcv, Evnt, Alloc>::_init_structures()
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
template <typename Tr, typename Vis, typename Subcv, typename Evnt,
          typename Alloc>
void Basic_sweep_line_2<Tr, Vis, Subcv, Evnt, Alloc>::_complete_sweep()
{
  CGAL_assertion(m_queue->empty());
  CGAL_assertion((m_statusLine.size() == 0));
  
  // Free all subcurve objects.
  for (unsigned int i = 0; i < m_num_of_subCurves; ++i)
    m_subCurveAlloc.destroy(m_subCurves + i);

  if (m_num_of_subCurves > 0)
    m_subCurveAlloc.deallocate(m_subCurves, m_num_of_subCurves);
}

//-----------------------------------------------------------------------------
// Initialize an event associated with a point.
//
template <typename Tr, typename Vis, typename Subcv, typename Evnt,
          typename Alloc>
void Basic_sweep_line_2<Tr, Vis, Subcv, Evnt, Alloc>::
_init_point(const Point_2& pt, Attribute type)
{
  // Create the event, or obtain an existing event in the queue.
  // Note that an isolated point does not have any boundary conditions.
  const std::pair<Event*, bool>& pair_res =
    _push_event(pt, type, ARR_INTERIOR, ARR_INTERIOR);

  bool is_new = pair_res.second;
  m_visitor->update_event(pair_res.first, pt, is_new);
}
  
//-----------------------------------------------------------------------------
// Initialize the events associated with an x-monotone curve.
//
template <typename Tr, typename Vis, typename Subcv, typename Evnt,
          typename Alloc>
void Basic_sweep_line_2<Tr, Vis, Subcv, Evnt, Alloc>::
_init_curve(const X_monotone_curve_2& curve, unsigned int index)
{
  // Construct an initialize a subcurve object.
  m_subCurveAlloc.construct(m_subCurves + index, m_masterSubcurve);
  
  (m_subCurves + index)->init(curve);
  
  // Create two events associated with the curve ends.
  _init_curve_end(curve, ARR_MAX_END, m_subCurves + index);
  _init_curve_end(curve, ARR_MIN_END, m_subCurves + index);
}

//-----------------------------------------------------------------------------
// Initialize an event associated with an x-monotone curve end.
//
template <typename Tr, typename Vis, typename Subcv, typename Evnt,
          typename Alloc>
void Basic_sweep_line_2<Tr, Vis, Subcv, Evnt, Alloc>::
_init_curve_end(const X_monotone_curve_2& cv, Arr_curve_end ind, Subcurve* sc)
{
  // Get the boundary conditions of the curve end.
  const Attribute  end_attr =
    (ind == ARR_MIN_END) ? Base_event::LEFT_END : Base_event::RIGHT_END;

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
      _push_event(pt, end_attr, ps_x, ps_y, sc) :
      _push_event(cv, ind, end_attr, ps_x, ps_y, sc);

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
template <typename Tr, typename Vis, typename Subcv, typename Evnt,
          typename Alloc>
void Basic_sweep_line_2<Tr, Vis, Subcv, Evnt, Alloc>::_handle_left_curves()
{ 
  CGAL_PRINT("Handling left curve" << std::endl;);
  
  m_is_event_on_above = false;
  
  if (! m_currentEvent->has_left_curves()) { 
    // The current event does not have any incident left subcurves, so we
    // should find a place for it in the status line (the function we call
    // update the m_status_line_insert_hint and m_is_event_on_above members).
    // We also notify the visitor on the new event we are about to handle.
    _handle_event_without_left_curves();

    if (m_currentEvent->is_closed()) {
      if (m_is_event_on_above) {
        // The current event is on the interior of existing curve on the
        // status line. Since the basic sweep does not allow intersections,
        // this is possible  only if the event is an isolated query point.
        CGAL_assertion(! m_currentEvent->has_right_curves() &&
                        m_currentEvent->is_query());
        
        //m_is_event_on_above = true;
        m_visitor->before_handle_event(m_currentEvent);
      }
      else
        m_visitor->before_handle_event(m_currentEvent);
    }
    else
      m_visitor->before_handle_event(m_currentEvent);
    
    // Nothing else to do (no left curves).
    return;
  }
  
  CGAL_PRINT("left curves before sorting: "<<"\n";);
  CGAL_SL_DEBUG(if (m_currentEvent->left_curves_begin() != 
                    m_currentEvent->left_curves_end() )
                {
                  m_currentEvent->Print();
                });

  // Use the status-line to sort all left subcurves incident to the current
  // event (no geometric comparisons are neede at all).
  _sort_left_curves();

  // Now the event is updated, with its left subcurved properly sorted, and
  // we can inform the visitor that we are about to handle this event.
  m_visitor->before_handle_event(m_currentEvent);
  
  CGAL_PRINT("left curves after sorting: " << std::endl);
  CGAL_SL_DEBUG(if (m_currentEvent->left_curves_begin() != 
                    m_currentEvent->left_curves_end())
                {
                  m_currentEvent->Print();
                });
  
  // Remove all left subcurves from the status line, and inform the visitor
  // that we are done handling these subcurves.
  Event_subcurve_iterator left_iter = m_currentEvent->left_curves_begin();
  
  while (left_iter != m_currentEvent->left_curves_end()) {
    Subcurve* left_sc = *left_iter;
 
    m_visitor->add_subcurve(left_sc->last_curve(), left_sc);
    ++left_iter;
    
    _remove_curve_from_status_line(left_sc);    
  }
  CGAL_PRINT( "Handling left curve END" << std::endl;);
}

//-----------------------------------------------------------------------------
// Handle an event that does not have any incident left curves.
//
template <typename Tr, typename Vis, typename Subcv, typename Evnt,
          typename Alloc>
void Basic_sweep_line_2<Tr, Vis, Subcv, Evnt, Alloc>::
_handle_event_without_left_curves()
{
  // Check if the event is a boundary event or not.
  const Arr_parameter_space  ps_x = m_currentEvent->parameter_space_in_x();
  const Arr_parameter_space  ps_y = m_currentEvent->parameter_space_in_y();

  if ((ps_x == ARR_INTERIOR) && (ps_y == ARR_INTERIOR)) {
    // The event is associated with a valid point - locate the position of
    // this point on the status line (note this point may be located on a
    // subcurve in the status line).
    const std::pair<Status_line_iterator, bool>& pair_res =
      m_statusLine.find_lower(m_currentEvent->point(), m_statusLineCurveLess);
    
    m_status_line_insert_hint = pair_res.first;
    m_is_event_on_above = pair_res.second;
    
    return;
  }
  
  // We have a boundary event, so we can easily locate a plave for it in the
  // status line.

  if (ps_x == ARR_LEFT_BOUNDARY) {
    // We are still sweeping the left boundary, so by the way we have ordered
    // the events in the queue, we know that the new event should be placed
    // above all other subcurves in the status line.
    m_status_line_insert_hint = m_statusLine.end();
  }
  else {
    // Note that an event with a positive boundary condition at x can only
    // represent a right end of a curve.
    CGAL_assertion(ps_x != ARR_RIGHT_BOUNDARY);

    // If the sign of the boundary in y is negative, the event should be
    // inserted below all other subcurves; if it is possitive, the event is
    // above all other subcurves.
    if (ps_y == ARR_BOTTOM_BOUNDARY) {
      m_status_line_insert_hint = m_statusLine.begin();
    }
    else {
      CGAL_assertion(ps_y == ARR_TOP_BOUNDARY);
      m_status_line_insert_hint = m_statusLine.end();
    }
  }
}

//-----------------------------------------------------------------------------
// Sort the left subcurves of an event point according to their order in
// their status line (no geometric comprasions are needed).
//
template <typename Tr, typename Vis, typename Subcv, typename Evnt,
          typename Alloc>
void Basic_sweep_line_2<Tr, Vis, Subcv, Evnt, Alloc>::_sort_left_curves()
{
  CGAL_assertion(m_currentEvent->has_left_curves());
  
  // Get the first curve associated with the event and its position on the
  // status line. We proceed from this position up the status line until
  // we encounter a subcurve that is not associated with the current event.
  Subcurve* curve = *(m_currentEvent->left_curves_begin());
  Status_line_iterator sl_iter = curve->hint();
  CGAL_assertion(*sl_iter == curve);
  // Look for the first curve in the vertical ordering that is also in the
  // left curve of the event  
  for (++sl_iter; sl_iter != m_statusLine.end(); ++sl_iter) {
    if (std::find(m_currentEvent->left_curves_begin(),
                  m_currentEvent->left_curves_end(), *sl_iter) ==
        m_currentEvent->left_curves_end())
      break;
  }
  Status_line_iterator end = sl_iter;
  
  sl_iter = curve->hint();
  if (sl_iter == m_statusLine.begin()) {
    // In case the lowest subcurve in the status line is associated with the
    // current event, we have the range of (sorted) subcurves ready. We
    // associate this range with the event, so the curves are now sorted
    // according to their vertical positions immediately to the left of the
    // event.
    m_currentEvent->replace_left_curves(sl_iter,end);
    return;
  }

  // Go down the status line until we encounter a subcurve that is not
  // associated with the current event.
  --sl_iter;
  for (;sl_iter != m_statusLine.begin(); --sl_iter) {
    if (std::find(m_currentEvent->left_curves_begin(),
                  m_currentEvent->left_curves_end(), *sl_iter) ==
        m_currentEvent->left_curves_end())
    {
      // Associate the sorted range of subcurves with the event.
      m_currentEvent->replace_left_curves(++sl_iter,end);
      return;
    }
  }

  // Check if the subcurve at the current iterator position should be
  // associated with the current event, and select the (sorted) range of
  // subcurves accordingly.
  if (std::find(m_currentEvent->left_curves_begin(),
                m_currentEvent->left_curves_end(), *sl_iter) ==
      m_currentEvent->left_curves_end())
    m_currentEvent->replace_left_curves(++sl_iter,end);
  else
    m_currentEvent->replace_left_curves(sl_iter,end);
}

//-----------------------------------------------------------------------------
// Handle the subcurves to the right of the current event point.
//
template <typename Tr, typename Vis, typename Subcv, typename Evnt,
          typename Alloc>
void Basic_sweep_line_2<Tr, Vis, Subcv, Evnt, Alloc>::_handle_right_curves()
{
  CGAL_PRINT("Handling right curves (");
  CGAL_SL_DEBUG(PrintEvent(m_currentEvent));
  CGAL_PRINT(")\n");
    
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
  Status_line_iterator sl_iter;

  while (curr != right_end) {
    CGAL_PRINT_INSERT(*curr);
    sl_iter = m_statusLine.insert_before(m_status_line_insert_hint, *curr);
    ((Subcurve*)(*curr))->set_hint(sl_iter);
    
    CGAL_SL_DEBUG(PrintStatusLine(););
    ++curr;
  }        
      
  CGAL_SL_DEBUG(PrintStatusLine());
}

//-----------------------------------------------------------------------------
// Add a subcurve to the right of an event point.
//
template <typename Tr, typename Vis, typename Subcv, typename Evnt,
          typename Alloc>
bool Basic_sweep_line_2<Tr, Vis, Subcv, Evnt, Alloc>::
_add_curve_to_right(Event* event, Subcurve* curve, bool /* overlap_exist */)
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
template <typename Tr, typename Vis, typename Subcv, typename Evnt,
          typename Alloc>
void Basic_sweep_line_2<Tr, Vis, Subcv, Evnt, Alloc>::
_remove_curve_from_status_line(Subcurve* sc)
{
  CGAL_PRINT("remove_curve_from_status_line\n";);
  CGAL_SL_DEBUG(PrintStatusLine(););
  CGAL_SL_DEBUG(sc->Print(););

  // Get the position of the subcurve on the status line.
  Status_line_iterator sl_iter = sc->hint(); 
  CGAL_assertion(sl_iter != m_statusLine.end());

  // The position of the next event can be right after the deleted subcurve.
  m_status_line_insert_hint = sl_iter; 
  ++m_status_line_insert_hint; 

  // Erase the subcurve from the status line.
  m_statusLine.erase(sl_iter);
  CGAL_PRINT("remove_curve_from_status_line Done\n";)
}

//-----------------------------------------------------------------------------
// Allocate an event object associated with a valid point.
//
template <typename Tr, typename Vis, typename Subcv, typename Evnt,
          typename Alloc>
typename Basic_sweep_line_2<Tr, Vis, Subcv, Evnt, Alloc>::Event*
Basic_sweep_line_2<Tr, Vis, Subcv, Evnt, Alloc>::
_allocate_event(const Point_2& pt, Attribute type,
                Arr_parameter_space ps_x, Arr_parameter_space ps_y)
{
  // Allocate the event.
  Event* e =  m_eventAlloc.allocate(1); 
  m_eventAlloc.construct(e, m_masterEvent);
  e->init(pt, type, ps_x, ps_y);

  // Insert it to the set of allocated events.
  m_allocated_events.insert(e);
  return e;
}

//-----------------------------------------------------------------------------
// Allocate an event at open boundary, 
// which is not associated with a valid point.
//
template <typename Tr, typename Vis, typename Subcv, typename Evnt,
          typename Alloc>
typename Basic_sweep_line_2<Tr, Vis, Subcv, Evnt, Alloc>::Event*
Basic_sweep_line_2<Tr, Vis, Subcv, Evnt, Alloc>::
_allocate_event_at_open_boundary(Attribute type,
                                 Arr_parameter_space ps_x, 
                                 Arr_parameter_space ps_y)
{
  Event* e =  m_eventAlloc.allocate(1); 
  m_eventAlloc.construct(e, m_masterEvent);
  e->init_at_open_boundary(type, ps_x, ps_y);

  m_allocated_events.insert(e);
  return e;
}

//-----------------------------------------------------------------------------
// Push a closed event point into the event queue.
//
template <typename Tr, typename Vis, typename Subcv, typename Evnt,
          typename Alloc>
std::pair<typename Basic_sweep_line_2<Tr, Vis, Subcv, Evnt, Alloc>::Event*,
          bool>
Basic_sweep_line_2<Tr, Vis, Subcv, Evnt, Alloc>::
_push_event(const Point_2& pt, Attribute type,
            Arr_parameter_space ps_x, Arr_parameter_space ps_y, Subcurve* sc)
{
  // Look for the point in the event queue.
  Event* e;  
  m_queueEventLess.set_parameter_space_in_x(ps_x);
  m_queueEventLess.set_parameter_space_in_y(ps_y);  
  
  const std::pair<Event_queue_iterator, bool>& pair_res =
    m_queue->find_lower(pt, m_queueEventLess);
  const bool exist = pair_res.second;
  if (! exist) {
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

  // If we are given a subcurve that the event represents one of its
  // endpoints, update the event and the subcurve records accordingly.
  // Note that this must be done before we actually insert the new event
  // into the event queue.
  if (sc != NULL) {
    if (type == Base_event::LEFT_END) {
      sc->set_left_event(e);
      _add_curve_to_right(e, sc);
    }
    else {
      CGAL_assertion(type == Base_event::RIGHT_END);
      sc->set_right_event(e);
      e->add_curve_to_left(sc);
    }
  }
 
  if (! exist) {
    // Insert the new event into the queue using the hint we got when we
    // looked for it.
    m_queue->insert_before(pair_res.first, e);
    CGAL_PRINT_NEW_EVENT(pt, e);
  }
  else {
    CGAL_PRINT_UPDATE_EVENT(pt, e);
  }
  
  // Return the resulting event and a flag indicating whether we have created
  // a new event.
  return (std::make_pair(e, !exist));
}

//-----------------------------------------------------------------------------
// Push an event point associated with a curve end into the event queue.
//
template <typename Tr, typename Vis, typename Subcv, typename Evnt,
          typename Alloc>
std::pair<typename Basic_sweep_line_2<Tr, Vis, Subcv, Evnt, Alloc>::Event*,
          bool>
Basic_sweep_line_2<Tr, Vis, Subcv, Evnt, Alloc>::
_push_event(const X_monotone_curve_2& cv, Arr_curve_end ind, Attribute type,
            Arr_parameter_space ps_x, Arr_parameter_space ps_y, Subcurve* sc)
{
  //cv has no member named 'base'
  //std::cout << "cv: " << cv.base() << std::endl;
  
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
    CGAL_assertion((e->parameter_space_in_x() == ps_x) &&
                   (e->parameter_space_in_y() == ps_y));

    e->set_attribute(type);
  }

  // If we are given a subcurve that the event represents one of its
  // endpoints, update the event and the subcurve records accordingly.
  // Note that this must be done before we actually insert the new event
  // into the event queue.
  if (sc != NULL) {
    if (type == Base_event::LEFT_END) {
      sc->set_left_event(e);
      _add_curve_to_right(e, sc);
    }
    else {
      CGAL_assertion(type == Base_event::RIGHT_END);
      sc->set_right_event(e);
      e->add_curve_to_left(sc);
    }
  }

  if (! exist) {
    // Insert the new event into the queue using the hint we got when we
    // looked for it.
    m_queue->insert_before(pair_res.first, e);
  }
  return (std::make_pair(e, !exist));
}

} //namespace CGAL

#endif
