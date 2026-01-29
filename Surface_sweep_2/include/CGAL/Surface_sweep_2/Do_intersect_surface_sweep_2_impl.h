// Copyright (c) 2025 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Efi Fogel        <efif@post.tau.ac.il>

#ifndef CGAL_SURFACE_SWEEP_2_DO_INTERSECT_SURFACE_SWEEP_2_IMPL_H
#define CGAL_SURFACE_SWEEP_2_IMPL_H

#include <CGAL/license/Surface_sweep_2.h>
#include <set>
#include <vector>

/*! \file
 *
 * Member-function definitions of the Do_intersect_surface_sweep_2 class-template.
 */

#include <algorithm>
#include <set>

namespace CGAL {
namespace Surface_sweep_2 {

/*! completes the sweep (complete the data structures).
 */
template <typename Vis>
void Do_intersect_surface_sweep_2<Vis>::_complete_sweep() {
  this->m_queue->clear();
  this->m_statusLine.clear();
  this->_clear();
}

/*! handles the subcurves to the left of the current event point.
 */
template <typename Vis>
void Do_intersect_surface_sweep_2<Vis>::_handle_left_curves() {
  this->m_is_event_on_above = false;

  if (! this->m_currentEvent->has_left_curves()) {
    // In case the current event has no left subcurves incident to it, we have
    // to locate a place for it in the status line.
    this->_handle_event_without_left_curves(Sides_category());

    if (this->m_is_event_on_above) {
      this->m_visitor->found_intersection();
      return;
    }

    // The event is not located on any subcurve.
    this->m_visitor->before_handle_event(this->m_currentEvent);
    return;
  }

  this->_sort_left_curves();
  this->m_visitor->before_handle_event(this->m_currentEvent);

  // Remove all left subcurves from the status line, and inform the visitor
  // that we are done handling these subcurves.
  Event_subcurve_iterator left_iter = this->m_currentEvent->left_curves_begin();
  while (left_iter != this->m_currentEvent->left_curves_end()) {
    Subcurve* left_sc = *left_iter;
    CGAL_assertion(left_sc->right_event() == this->m_currentEvent);
    this->m_visitor->add_subcurve(left_sc->last_curve(), left_sc);
    ++left_iter;
    if (_remove_curve_from_status_line(left_sc)) return;
  }
}

/*! handles the subcurves to the right of the current event point.
 */
template <typename Vis>
void Do_intersect_surface_sweep_2<Vis>::_handle_right_curves() {
  for (auto it = this->m_currentEvent->right_curves_begin(); it != this->m_currentEvent->right_curves_end(); ++it) {
    Subcurve* subcurve = *it;
    subcurve->reset_left_event();
  }

  if (! this->m_currentEvent->has_right_curves()) return;

  auto curr = this->m_currentEvent->right_curves_begin();
  auto right_end = this->m_currentEvent->right_curves_end();

  auto sl_it = this->m_statusLine.insert_before(this->m_status_line_insert_hint, *curr);
  Subcurve* sc = *curr;
  sc->set_hint(sl_it);

  if (sl_it != this->m_statusLine.begin()) {
    auto prev = sl_it;
    --prev;
    if (_do_intersect(*prev, *sl_it)) return;
  }

  auto prev = curr;
  ++curr;
  while (curr != right_end) {
    sl_it = this->m_statusLine.insert_before(this->m_status_line_insert_hint, *curr);

    Subcurve* sc = *curr;
    sc->set_hint(sl_it);

    // If the two curves used to be neighbors before, we do not need to
    // intersect them again.
    if (! this->m_currentEvent->are_left_neighbors(*curr, *prev)) {
      if (_do_intersect(*prev, *curr)) return;
    }

    prev = curr;
    ++curr;
  }

  //the next Subcurve at the status line
  ++sl_it;
  if (sl_it != this->m_statusLine.end()) {
    if (_do_intersect(*prev, *sl_it)) return;
  }
}

/*! adds a subcurve to the right of an event point.
 */
template <typename Vis>
bool Do_intersect_surface_sweep_2<Vis>::_add_curve_to_right(Event* event, Subcurve* curve) {
  std::pair<bool, Event_subcurve_iterator> pair_res = event->add_curve_to_right(curve, this->m_traits);
  if (! pair_res.first) return false;    // no overlap
  this->m_visitor->found_intersection();
  return true;
}

/*! removes a curve from the status line.
 */
template <typename Vis>
bool Do_intersect_surface_sweep_2<Vis>::_remove_curve_from_status_line(Subcurve* sc) {
  auto sl_it = sc->hint();
  this->m_status_line_insert_hint = sl_it;
  ++(this->m_status_line_insert_hint);
  sc->set_hint(this->m_statusLine.end());

  // The subcurve is removed for good from the status line. We need
  // to check for intersection between its two neighbors (below and above)
  CGAL_assertion(sl_it != this->m_statusLine.end());
  auto last_it = this->m_statusLine.end();
  --last_it;

  if ((sl_it != this->m_statusLine.begin()) && (sl_it != last_it)) {
    auto prev = sl_it;
    --prev;
    auto next = sl_it;
    ++next;
    if (_do_intersect(*prev, *next)) return true;
  }
  this->m_statusLine.erase(sl_it);
  return false;
}

/*! computes intersections between the two given curves.
 */
template <typename Vis>
bool Do_intersect_surface_sweep_2<Vis>::_do_intersect(Subcurve* c1, Subcurve* c2) {
  CGAL_assertion(c1 != c2);
  auto do_intersect = this->m_traits->do_intersect_2_object();
  if (! do_intersect(c1->last_curve(), c2->last_curve(), m_closed)) return false;
  this->m_visitor->found_intersection();
  return true;
}

} // namespace Surface_sweep_2
} // namespace CGAL

#endif
