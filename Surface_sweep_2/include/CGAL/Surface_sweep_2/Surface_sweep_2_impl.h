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
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s) : Baruch Zukerman <baruchzu@post.tau.ac.il>
//             Efi Fogel       <efif@post.tau.ac.il>
//               (based on old version by Tali Zvi)

#ifndef CGAL_SURFACE_SWEEP_2_IMPL_H
#define CGAL_SURFACE_SWEEP_2_IMPL_H

#include <CGAL/license/Surface_sweep_2.h>

/*! \file
 *
 * Member-function definitions of the Surface_sweep_2 class-template.
 */

namespace CGAL {
namespace Surface_sweep_2 {

//-----------------------------------------------------------------------------
// Initialize the data structures for the sweep-line algorithm.
//
template <typename Vis>
void Surface_sweep_2<Vis>::_init_structures()
{
  // Initailize the structures maintained by the base sweep-line class.
  Base::_init_structures();
}

//-----------------------------------------------------------------------------
// Complete the sweep (complete the data structures).
//
template <typename Vis>
void Surface_sweep_2<Vis>::_complete_sweep()
{
  CGAL_SS_PRINT_START_EOL("completing the sweep");

  // Complete the sweep process using base sweep-line class.
  Base::_complete_sweep();

  // Clean the set of curve pairs for which we have computed intersections.
  m_curves_pair_set.clear();

  // Free all overlapping subcurves we have created.
  Subcurve_iterator   itr;
  for (itr = m_overlap_subCurves.begin(); itr != m_overlap_subCurves.end();
       ++itr)
  {
#ifdef CGAL_CXX11
    std::allocator_traits<Subcurve_alloc>::destroy(this->m_subCurveAlloc, *itr);
#else
    this->m_subCurveAlloc.destroy(*itr);
#endif
    this->m_subCurveAlloc.deallocate(*itr, 1);
  }

  m_overlap_subCurves.clear();

  CGAL_SS_PRINT_END_EOL("completing the sweep");
}

//-----------------------------------------------------------------------------
// Handle the subcurves to the left of the current event point.
//
template <typename Vis>
void Surface_sweep_2<Vis>::_handle_left_curves()
{
  CGAL_SS_PRINT_START("handling left curves at (");
  CGAL_SS_DEBUG(this->PrintEvent(this->m_currentEvent));
  CGAL_SS_PRINT_TEXT(")");
  CGAL_SS_PRINT_EOL();

  this->m_is_event_on_above = false;

  if (! this->m_currentEvent->has_left_curves()) {
    // In case the current event has no left subcurves incident to it, we have
    // to locate a place for it in the status line.
    CGAL_SS_PRINT_TEXT("Handling case: no left curves");
    CGAL_SS_PRINT_EOL();
    this->_handle_event_without_left_curves();

    // Status_line_iterator sl_pos = this->m_status_line_insert_hint;

    if (this->m_is_event_on_above) {
      CGAL_SS_PRINT_TEXT("The event is on a curve in the status line");
      CGAL_SS_PRINT_EOL();


      // Obtain the subcurve that contains the current event
      Subcurve* sc = *(this->m_status_line_insert_hint);

      // The current event point starts at the interior of a subcurve that
      // already exists in the status line (this may also indicate an overlap).
      if (! this->m_currentEvent->has_right_curves()) {
        // The event is an isolated point.
        if (this->m_currentEvent->is_query()) {
          // In case of a query point, just notify the visitor about it.
          this->m_is_event_on_above = true;
          this->m_visitor->before_handle_event(this->m_currentEvent);
          return;
        }

        // In case of an isolated action point, mark the point at a "weak"
        // intersection.
        CGAL_assertion(this->m_currentEvent->is_action());
        this->m_currentEvent->set_weak_intersection();
        this->m_visitor->update_event(this->m_currentEvent, sc);

        this->m_currentEvent->add_curve_to_left(sc);
        this->m_currentEvent->push_back_curve_to_right(sc);
      }
      else
      {
        this->m_currentEvent->push_back_curve_to_left(sc);
        this->m_currentEvent->set_weak_intersection();
        this->m_visitor->update_event(this->m_currentEvent, sc);
        _add_curve_to_right(this->m_currentEvent, sc);
      }

      // sc is now on the left
      CGAL_SS_PRINT_TEXT("Event after update:");
      CGAL_SS_PRINT_EOL();
      CGAL_SS_PRINT_EVENT_INFO(this->m_currentEvent);
      CGAL_SS_PRINT_EOL();
      CGAL_assertion(std::distance(this->m_currentEvent->left_curves_begin(),
                                   this->m_currentEvent->left_curves_end())==1);
    }
    else {
      // The event is not located on any subcurve.
      this->m_visitor->before_handle_event(this->m_currentEvent);
      CGAL_SS_PRINT_END_EOL("handling left curves");
      return;
    }
  }

  CGAL_SS_PRINT_TEXT("left curves before sorting:");
  CGAL_SS_PRINT_EOL();
  CGAL_SS_DEBUG(if (this->m_currentEvent->left_curves_begin() !=
                    this->m_currentEvent->left_curves_end())
                { this->print_event_info(this->m_currentEvent); });
  _fix_overlap_subcurves();
  this->_sort_left_curves();
  this->m_visitor->before_handle_event(this->m_currentEvent);

  CGAL_SS_PRINT_TEXT("left curves after sorting:");
  CGAL_SS_PRINT_EOL();
  CGAL_SS_DEBUG(if (this->m_currentEvent->left_curves_begin() !=
                    this->m_currentEvent->left_curves_end() )
                { this->print_event_info(this->m_currentEvent); });

  // Check if the curve should be removed for good.
  bool remove_for_good = false;

  Event_subcurve_iterator left_iter =
    this->m_currentEvent->left_curves_begin();
  while (left_iter != this->m_currentEvent->left_curves_end()) {
    Subcurve* leftCurve = *left_iter;

    if ((Event*)leftCurve->right_event() == this->m_currentEvent) {
      // we are done with that subcurve (current event point is his right
      // end point) so we remove it from the status line for good.
      remove_for_good = true;
      this->m_visitor->add_subcurve(leftCurve->last_curve(), leftCurve);
    }
    else {
      // curren event splits the subcurve.
      const X_monotone_curve_2& lastCurve = leftCurve->last_curve();
      this->m_traits->split_2_object()(lastCurve, this->m_currentEvent->point(),
                                       sub_cv1, sub_cv2);
      this->m_visitor->add_subcurve(sub_cv1, leftCurve);
      leftCurve->set_last_curve(sub_cv2);
    }
    ++left_iter;

    //remove curve from the status line (also checks intersection
    //between the neighbouring curves,only if the curve is removed for good)
    _remove_curve_from_status_line(leftCurve, remove_for_good);
  }

  CGAL_SS_PRINT_END_EOL("handling left curves");
}

//-----------------------------------------------------------------------------
// Handle the subcurves to the right of the current event point.
//
template <typename Vis>
void Surface_sweep_2<Vis>::_handle_right_curves()
{
  CGAL_SS_PRINT_START("handling right curves at (");
  CGAL_SS_DEBUG(this->PrintEvent(this->m_currentEvent));
  CGAL_SS_PRINT_TEXT(")");
  CGAL_SS_PRINT_EOL();

  if (! this->m_currentEvent->has_right_curves()) {
    CGAL_SS_PRINT_END_EOL("handling right curves");
    return;
  }

  // Some overlapping curves that are on the right might not be
  // on the left (in case of overlap with other curves having common
  // ancesters). Consequently, an overlapping curve on the right might
  // not have been split.
  // SL_SAYS: we should be able to do something better without geometric test
  // (like having an additional boolean in subcurve that we can set in
  //  _handle_left_curves() after a split)
  for(  Event_subcurve_iterator cit = this->m_currentEvent->right_curves_begin(),
                                cit_end = this->m_currentEvent->right_curves_end();
                                cit!=cit_end; ++cit)
  {
    if ( (*cit)->originating_subcurve1()!=NULL &&
         (Event*)(*cit)->left_event()!=this->m_currentEvent )
    {
      // split the subcurve.
      const X_monotone_curve_2& lastCurve = (*cit)->last_curve();

      if (!this->m_traits->equal_2_object()(
              this->m_traits->construct_min_vertex_2_object()(lastCurve),
              this->m_currentEvent->point()))
      {
        this->m_traits->split_2_object()(lastCurve, this->m_currentEvent->point(),
                                       sub_cv1, sub_cv2);
        (*cit)->set_last_curve(sub_cv2);
      }
    }
  }

  // Loop over the curves to the right of the status line and handle them:
  // - If we are at the beginning of the curve, we insert it to the status
  //   line, then we look if it intersects any of its neighbors.
  // - If we are at an intersection point between two curves, we add them
  //   to the status line and attempt to intersect them with their neighbors
  // - We also check to see if the two intersect again to the right of the
  //   point.

  Event_subcurve_iterator currentOne =
    this->m_currentEvent->right_curves_begin();
  Event_subcurve_iterator rightCurveEnd =
    this->m_currentEvent->right_curves_end();

  CGAL_SS_PRINT_INSERT(*currentOne);
  Status_line_iterator slIter =
    this->m_statusLine.insert_before(this->m_status_line_insert_hint,
                                     *currentOne);
  Subcurve* sc = *currentOne;
  sc->set_hint(slIter);

  CGAL_SS_PRINT_STATUS_LINE();
  if (slIter != this->m_statusLine.begin()) {
    //  get the previous curve in the y-str
    Status_line_iterator prev = slIter; --prev;
    _intersect(*prev, *slIter);
  }

  Event_subcurve_iterator prevOne = currentOne;
  ++currentOne;
  while (currentOne != rightCurveEnd) {
    CGAL_SS_PRINT_INSERT(*currentOne);
    slIter = this->m_statusLine.insert_before(this->m_status_line_insert_hint,
                                              *currentOne);

    Subcurve* sc = *currentOne;
    sc->set_hint(slIter);

    CGAL_SS_PRINT_STATUS_LINE();

    // If the two curves used to be neighbours before, we do not need to
    // intersect them again.
    if (!this->m_currentEvent->are_left_neighbours(*currentOne, *prevOne))
      _intersect(*prevOne, *currentOne);

    prevOne = currentOne;
    ++currentOne;
  }

  CGAL_SS_PRINT_STATUS_LINE();

  //the next Subcurve at the status line
  ++slIter;
  if (slIter != this->m_statusLine.end()) _intersect(*prevOne, *slIter);

  CGAL_SS_PRINT_END_EOL("handling right curves");
}

//-----------------------------------------------------------------------------
// Add a subcurve to the right of an event point.
//
template <typename Vis>
bool Surface_sweep_2<Vis>::_add_curve_to_right(Event* event, Subcurve* curve)
{
  CGAL_SS_PRINT_START("adding a Curve to the right of (");
  CGAL_SS_DEBUG(this->PrintEvent(event));
  CGAL_SS_PRINT_TEXT(") ");
  CGAL_SS_PRINT_CURVE(curve);
  CGAL_SS_PRINT_EOL();

  Event_subcurve_iterator iter;
  for (iter = event->right_curves_begin(); iter != event->right_curves_end();
       ++iter)
  {
    CGAL_SS_PRINT_CURVE(*iter);
    CGAL_SS_PRINT_EOL();
    if ((*iter)->are_all_leaves_contained(curve)) {
      CGAL_SS_PRINT_END_EOL("adding a Curve to the right (curve exists)");
      return false;
    }

    if (curve->are_all_leaves_contained(*iter)) {
      *iter = curve;    // replace the current curve with the new one.
      CGAL_SS_PRINT_END_EOL
        ("adding a Curve to the right (curve partially overlaps)");
      return false;
    }

    CGAL_assertion(!(curve)->has_same_leaves(*iter));
  }
  std::pair<bool, Event_subcurve_iterator> pair_res =
    event->add_curve_to_right(curve, this->m_traits);

  if (! pair_res.first) {
    // No overlap occurs.
    CGAL_SS_PRINT_END_EOL("adding a Curve to the right (no overlap)");
    return false;
  }

  // a new overlap needs to be computed
  _intersect(static_cast<Subcurve*>(curve),
             static_cast<Subcurve*>(*(pair_res.second)),
             event);
  // SL_SAYS: we need to be sure that the overlapping curve generated by
  //          `curve` and `*(pair_res.second)` has not yet been computed or
  //          that the overlapping curve will be added to the right of
  //          `event` because of an adjacency in the status line.

  // Indicate that an overlap has occured:
  CGAL_SS_PRINT_END_EOL("adding a Curve to the right (overlap)");
  return true;
}

//-----------------------------------------------------------------------------
// Remove a curve from the status line.
//
template <typename Vis>
void Surface_sweep_2<Vis>::_remove_curve_from_status_line(Subcurve* leftCurve,
                                                          bool remove_for_good)
{
  CGAL_SS_PRINT_START("removing a curve from the status line, ");
  CGAL_SS_PRINT_CURVE(leftCurve);
  CGAL_SS_PRINT_EOL();
  CGAL_SS_PRINT_STATUS_LINE();

  Status_line_iterator sliter = leftCurve->hint();
  this->m_status_line_insert_hint = sliter;
  ++(this->m_status_line_insert_hint);
  leftCurve->set_hint(this->m_statusLine.end());

  if (! remove_for_good) {
    // the subcurve is not removed for good, so we dont need to intersect
    // its neighbours after its removal.
    CGAL_SS_PRINT_ERASE(*sliter);
    this->m_statusLine.erase(sliter);
    CGAL_SS_PRINT_END_EOL("Removing a curve from the status line");
    return;
  }

  // the subcurve will be removed for good from the stauts line, we need
  // to check for intersection between his two neighbours (below and above him)
  // but we need to make sure that its not the first or last subcurve
  // at the status line.
  CGAL_assertion(sliter != this->m_statusLine.end());
  Status_line_iterator lastOne = this->m_statusLine.end();
  --lastOne;

  if (sliter != this->m_statusLine.begin() && sliter != lastOne) {
    Status_line_iterator prev = sliter; --prev;
    Status_line_iterator next = sliter; ++next;

    // intersect *next with  *prev
    _intersect(*prev, *next);
  }
  CGAL_SS_PRINT_ERASE(*sliter);
  this->m_statusLine.erase(sliter);
  CGAL_SS_PRINT_END_EOL("removing a curve from the status line");
}

//-----------------------------------------------------------------------------
// Compute intersections between the two given curves.
//
template <typename Vis>
  void Surface_sweep_2<Vis>::_intersect(Subcurve* c1, Subcurve* c2,
                                        Event* event_for_overlap)
{
  CGAL_SS_PRINT_START("computing intersection of ");
  CGAL_SS_PRINT_CURVE(c1);
  CGAL_SS_PRINT_TEXT(" and ");
  CGAL_SS_PRINT_CURVE(c2);
  CGAL_SS_PRINT_EOL();

  typedef typename Geometry_traits_2::Multiplicity      Multiplicity;

  CGAL_assertion(c1 != c2);

  // look up for (c1,c2) in the table and insert if doesnt exist
  Curve_pair cv_pair(c1,c2);
  if (! (m_curves_pair_set.insert(cv_pair)).second) {
    CGAL_SS_PRINT_END_EOL("computing intersection (already computed)");
    return;  //the curves have already been checked for intersection
  }

  // handle overlapping curves with common ancesters
  Subcurve_vector all_leaves_diff;
  Subcurve* first_parent=NULL;
  if (c1->originating_subcurve1()!=NULL || c2->originating_subcurve2()!=NULL)
  {
    // get the subcurve leaves of c1 and of c2. Then extract from the smallest set
    // the subcurves leaves that are not in the other one. If empty, it means that
    // a subcurves is completely contained in another one.
    first_parent = c1;
    Subcurve* second_parent = c2;

    Subcurve_vector all_leaves_first;
    Subcurve_vector all_leaves_second;
    first_parent->all_leaves(std::back_inserter(all_leaves_first));
    second_parent->all_leaves(std::back_inserter(all_leaves_second));
    if (all_leaves_second.size() > all_leaves_first.size())
    {
      std::swap(first_parent,second_parent);
      std::swap(all_leaves_first,all_leaves_second);
    }

    CGAL_assertion(!all_leaves_first.empty() && !all_leaves_second.empty());

    std::sort(all_leaves_first.begin(), all_leaves_first.end());
    std::sort(all_leaves_second.begin(), all_leaves_second.end());

    std::set_difference(all_leaves_second.begin(), all_leaves_second.end(),
                        all_leaves_first.begin(), all_leaves_first.end(),
                        std::back_inserter(all_leaves_diff));

    if (all_leaves_second.size()==all_leaves_diff.size())
      all_leaves_diff.clear(); // clear so that it is not used by _create_overlapping_curve()
    else
      if (all_leaves_diff.empty())
      {
        CGAL_SS_PRINT_TEXT("One overlapping curve entirely contains the other one");
        CGAL_SS_PRINT_EOL();

        Event* left_event = (Event*) first_parent->left_event();
        Event* right_event = (Event*) first_parent->right_event();

        if (!second_parent->is_start_point(left_event))
          left_event->add_curve_to_left(second_parent);
        else
          left_event->remove_curve_from_right(second_parent);

        CGAL_SS_PRINT_CURVE(c1);
        CGAL_SS_PRINT_TEXT(" + ");
        CGAL_SS_PRINT_CURVE(c2);
        CGAL_SS_PRINT_TEXT(" => ");
        CGAL_SS_PRINT_EOL();
        CGAL_SS_PRINT_TEXT("  ");
        CGAL_SS_PRINT_CURVE(first_parent);
        CGAL_SS_PRINT_EOL();

        // Remove second_parent from the left curves of the right end
        // and add it on the right otherwise
        if (second_parent->is_end_point(right_event))
          right_event->remove_curve_from_left(second_parent);
        else
          _add_curve_to_right(right_event, second_parent);

        // add the overlapping curve kept of the right of the left end
        _add_curve_to_right(left_event, first_parent);
        right_event->add_curve_to_left(first_parent);

        this->m_visitor->found_overlap(c1, c2, first_parent);

        CGAL_SS_PRINT_END_EOL("computing intersection");
        return;
      }
      else{
        CGAL_SS_PRINT_TEXT("Overlap with common ancestors");
        CGAL_SS_PRINT_EOL();

        // iteratively create the final overlapping (geometric) curve.
        // This is needed rather than simply computing the intersection of
        // the last curves of first_parent and second_parent as some traits
        // classes (such as Arr_curve_data_traits_2) override the Intersect_2
        // functor and expects the curve to have no common ancesters
        // (Arr_curve_data_traits_2 is used in the testsuite to sum up
        //  the overlapping degree of a curve)
        X_monotone_curve_2 xc = first_parent->last_curve();
        for (typename Subcurve_vector::iterator sc_it=all_leaves_diff.begin();
             sc_it!=all_leaves_diff.end(); ++sc_it)
        {
          std::vector<CGAL::Object> inter_res;
          ;

          this->m_traits->intersect_2_object()(xc,
                                               (*sc_it)->last_curve(),
                                               vector_inserter(inter_res));
          CGAL_assertion(inter_res.size()==1);
          CGAL_assertion( CGAL::object_cast< X_monotone_curve_2 >(&inter_res.front())!=NULL );
          xc = *CGAL::object_cast< X_monotone_curve_2 >(&inter_res.front());
        }
        _create_overlapping_curve(xc, c1 , c2, all_leaves_diff, first_parent, event_for_overlap);
        CGAL_SS_PRINT_END_EOL("computing intersection (overlap with common ancestors)");
        return;
      }
  }

  // do compute the intersection of the two curves
  vector_inserter vi(m_x_objects) ;
  vector_inserter vi_end(m_x_objects);

  vi_end =
    this->m_traits->intersect_2_object()(c1->last_curve(), c2->last_curve(), vi);

  if (vi == vi_end) {
    CGAL_SS_PRINT_END_EOL("Computing intersection (no intersection)");
    return; // no intersection at all
  }

  // The two subCurves may start at the same point, in that case we ignore the
  // first intersection point.

  const Arr_parameter_space ps_x1 =
    this->m_traits->parameter_space_in_x_2_object()(c1->last_curve(),
                                                    ARR_MIN_END);
  const Arr_parameter_space ps_y1 =
    this->m_traits->parameter_space_in_y_2_object()(c1->last_curve(),
                                                    ARR_MIN_END);
  const Arr_parameter_space ps_x2 =
    this->m_traits->parameter_space_in_x_2_object()(c2->last_curve(),
                                                    ARR_MIN_END);
  const Arr_parameter_space ps_y2 =
    this->m_traits->parameter_space_in_y_2_object()(c2->last_curve(),
                                                    ARR_MIN_END);

  if ((ps_x1 == ps_x2) && (ps_y1 == ps_y2) &&
      ((ps_x1 != ARR_INTERIOR) || (ps_y1 != ARR_INTERIOR)) &&
      this->m_traits->is_closed_2_object()(c1->last_curve(), ARR_MIN_END) &&
      this->m_traits->is_closed_2_object()(c2->last_curve(), ARR_MIN_END))
  {
    if ( object_cast<std::pair<Point_2, Multiplicity> >(&(*vi)) != NULL
         && this->m_traits->equal_2_object()
        (this->m_traits->construct_min_vertex_2_object()(c1->last_curve()),
         this->m_traits->construct_min_vertex_2_object()(c2->last_curve())))
    {
      CGAL_SS_PRINT_TEXT("Skipping common left endpoint on boundary ...");
      CGAL_SS_PRINT_EOL();
      ++vi;
    }
  }

  // If the two subcurves have a common right-event, and the last intersection
  // object is a point, we can ignore last intersection (note that in case of
  // an overlap that ends at the common endpoint, we definately want to keep
  // the intersection object).
  if (c1->right_event() == c2->right_event()) {
    vector_inserter vi_last = vi_end;

    --vi_last;
    if (object_cast<std::pair<Point_2, Multiplicity> >(&(*vi_last)) != NULL) {
      CGAL_SS_PRINT_TEXT("Skipping common right endpoint...");
      CGAL_SS_PRINT_EOL();
      --vi_end;
    }
  }
  else {
    // In case both right curve-ends have boundary conditions and are not
    // open, check whether the right endpoints are the same. If they are,
    // skip the last intersection point.
    const Arr_parameter_space ps_x1 =
      this->m_traits->parameter_space_in_x_2_object()(c1->last_curve(),
                                                      ARR_MAX_END);
    const Arr_parameter_space ps_y1 =
      this->m_traits->parameter_space_in_y_2_object()(c1->last_curve(),
                                                      ARR_MAX_END);
    const Arr_parameter_space ps_x2 =
      this->m_traits->parameter_space_in_x_2_object()(c2->last_curve(),
                                                      ARR_MAX_END);
    const Arr_parameter_space ps_y2 =
      this->m_traits->parameter_space_in_y_2_object()(c2->last_curve(),
                                                      ARR_MAX_END);

    if ((ps_x1 == ps_x2) && (ps_y1 == ps_y2) &&
        ((ps_x1 != ARR_INTERIOR) || (ps_y2 != ARR_INTERIOR)) &&
        this->m_traits->is_closed_2_object()(c1->last_curve(), ARR_MAX_END) &&
        this->m_traits->is_closed_2_object()(c2->last_curve(), ARR_MAX_END))
    {
      if (this->m_traits->equal_2_object()
          (this->m_traits->construct_max_vertex_2_object()(c1->last_curve()),
           this->m_traits->construct_max_vertex_2_object()(c2->last_curve())))
      {
        vector_inserter vi_last = vi_end;

        --vi_last;
        if (object_cast<std::pair<Point_2, Multiplicity> >(&(*vi_last)) != NULL)
        {
          CGAL_SS_PRINT_TEXT("Skipping common right endpoint on boundary...");
          CGAL_SS_PRINT_EOL();
          --vi_end;
        }
      }
    }
  }

  const std::pair<Point_2, Multiplicity>* xp_point;

  // Efi: why not skipping in a loop?check only one (that is, why not in a loop)?
  if (vi != vi_end) {
    xp_point = object_cast<std::pair<Point_2, Multiplicity> >(&(*vi));
    if (xp_point != NULL) {
      // Skip the intersection point if it is not larger than the current
      // event.
      if (this->m_queueEventLess(xp_point->first, this->m_currentEvent) !=
          LARGER)
      {
        ++vi;
      }
    }
  }

  for (; vi != vi_end; ++vi) {
    const X_monotone_curve_2* icv;
    unsigned int multiplicity = 0;

    xp_point = object_cast<std::pair<Point_2, Multiplicity> >(&(*vi));
    if (xp_point != NULL) {
      Point_2 xp = xp_point->first;
      multiplicity = xp_point->second;
      CGAL_SS_PRINT_TEXT("Found an intersection point");
      CGAL_SS_PRINT_EOL();
      _create_intersection_point(xp, multiplicity, c1, c2);
    }
    else {
      icv = object_cast<X_monotone_curve_2>(&(*vi));
      CGAL_assertion(icv != NULL);
      CGAL_SS_PRINT_TEXT("Found an overlap");
      CGAL_SS_PRINT_EOL();
      _create_overlapping_curve(*icv, c1 , c2, all_leaves_diff, first_parent, event_for_overlap);
    }
  }

  CGAL_SS_PRINT_END_EOL("computing intersection");
}

//-----------------------------------------------------------------------------
// Create an intersection-point event between two curves.
//
template <typename Vis>
void Surface_sweep_2<Vis>::_create_intersection_point(const Point_2& xp,
                                                      unsigned int multiplicity,
                                                      Subcurve*& c1,
                                                      Subcurve*& c2)
{
  CGAL_SS_PRINT_START_EOL("creating an intersection point between");
  CGAL_SS_PRINT_CURVE(c1);
  CGAL_SS_PRINT_EOL();
  CGAL_SS_PRINT_CURVE(c2);
  CGAL_SS_PRINT_EOL();

  // insert the event and check if an event at this point already exists.
  const std::pair<Event*, bool>& pair_res =
    this->_push_event(xp, Event::DEFAULT, ARR_INTERIOR, ARR_INTERIOR);

  Event* e = pair_res.first;
  if (pair_res.second) {
    // a new event is created , which indicates that the intersection point
    // cannot be one of the end-points of two curves
    CGAL_SS_PRINT_TEXT("A new event is created .. (");
    CGAL_SS_PRINT(xp);
    CGAL_SS_PRINT_TEXT(")");
    CGAL_SS_PRINT_EOL();

    e->set_intersection();

    this->m_visitor->update_event(e, c1, c2, true);
    e->push_back_curve_to_left(c1);
    e->push_back_curve_to_left(c2);

    // Act according to the multiplicity:
    if (multiplicity == 0) {
      // The multiplicity of the intersection point is unkown or undefined:
      _add_curve_to_right(e, c1);
      _add_curve_to_right(e, c2);
      if (e->is_right_curve_bigger(c1, c2, this->m_traits)) std::swap(c1, c2);
    }
    else {
      if ((multiplicity % 2) == 1) {
        // The mutiplicity of the intersection point is odd: Swap their
        // order to the right of this point.
        std::swap(c1,c2);
        e->add_curve_pair_to_right(c1, c2);
      }
      else {
        // The mutiplicity of the intersection point is even, so they
        // maintain their order to the right of this point.
        CGAL_assertion((multiplicity % 2) == 0);
        e->add_curve_pair_to_right(c1, c2);
      }
    }
  }
  else {
    // The event already exists, so we need to update it accordingly
    CGAL_SS_PRINT_TEXT("Event already exists, updating.. (");
    CGAL_SS_PRINT(xp);
    CGAL_SS_PRINT_TEXT(")");
    CGAL_SS_PRINT_EOL();

    if (!c1->is_start_point(e)) e->add_curve_to_left(c1);
    if (!c2->is_start_point(e)) e->add_curve_to_left(c2);

    if (!c1->is_end_point(e) && !c2->is_end_point(e)) {
      _add_curve_to_right(e, c1);
      _add_curve_to_right(e, c2);
      e->set_intersection();
      this->m_visitor->update_event(e, c1, c2, false);

      if (multiplicity==0)
      {
        if (e->is_right_curve_bigger(c1, c2, this->m_traits)) std::swap(c1, c2);
      }
      else
      {
        if (multiplicity%2==1) std::swap(c1, c2);
      }
    }
    else {
      if (!c1->is_end_point(e) && c2->is_end_point(e)) {
        _add_curve_to_right(e, c1);
        e->set_weak_intersection();
        this->m_visitor->update_event(e, c1);
      }
      else {
        if (c1->is_end_point(e) && !c2->is_end_point(e)) {
          _add_curve_to_right(e, c2);
          e->set_weak_intersection();
          this->m_visitor->update_event(e, c2);
        }
      }
    }


    CGAL_SS_PRINT_EVENT_INFO(e);
  }

  CGAL_SS_PRINT_END_EOL("Creating an intersection point");
}

//-----------------------------------------------------------------------------
// Fix overlap Subcurves before handling the current event.
//
template <typename Vis>
void Surface_sweep_2<Vis>::_fix_overlap_subcurves()
{
  CGAL_SS_PRINT_START_EOL("fixing overlap subcurves");

  CGAL_assertion(this->m_currentEvent->has_left_curves());

  Event_subcurve_iterator iter = this->m_currentEvent->left_curves_begin();

  //special treatment for Subcurves that store overlaps
  while (iter != this->m_currentEvent->left_curves_end()) {
    Subcurve* leftCurve = *iter;

    // we check if the subcurve store overlap and current event is its
    // right end point.
    if ((Event*)leftCurve->right_event() == this->m_currentEvent) {
      if (leftCurve->originating_subcurve1() != NULL) {
        Subcurve* orig_sc_1 = leftCurve->originating_subcurve1();
        Subcurve* orig_sc_2 = leftCurve->originating_subcurve2();

        _fix_finished_overlap_subcurve(orig_sc_1);
        _fix_finished_overlap_subcurve(orig_sc_2);
      }
    }
    ++iter;
  }

  CGAL_SS_PRINT_END_EOL("Fixing overlap subcurves");
}

template <typename Vis>
void Surface_sweep_2<Vis>::
_create_overlapping_curve(const X_monotone_curve_2& overlap_cv,
                          Subcurve*& c1 , Subcurve*& c2,
                          const Subcurve_vector& all_leaves_diff,
                          Subcurve* first_parent,
                          Event* event_on_overlap)
{
  // An overlap occurs:
  CGAL_SS_PRINT_START_EOL("creating an overlapping curve");

  // Get the left end of overlap_cv.
  Event* left_event;
  Arr_parameter_space  ps_x_l =
    this->m_traits->parameter_space_in_x_2_object()(overlap_cv, ARR_MIN_END);
  Arr_parameter_space  ps_y_l =
    this->m_traits->parameter_space_in_y_2_object()(overlap_cv, ARR_MIN_END);
  if ((ps_x_l != ARR_INTERIOR) || (ps_y_l != ARR_INTERIOR)) {
    CGAL_assertion(c1->left_event() == c2->left_event());
    left_event=(Event*)(c1->left_event());
  }
  else{
    Point_2 left_end = this->m_traits->construct_min_vertex_2_object()(overlap_cv);
    left_event = this->_push_event(left_end, Event::DEFAULT, ARR_INTERIOR, ARR_INTERIOR).first;
  }

  // Get the right end of overlap_cv.
  Event* right_event;
  Arr_parameter_space  ps_x_r =
    this->m_traits->parameter_space_in_x_2_object()(overlap_cv, ARR_MAX_END);
  Arr_parameter_space ps_y_r =
    this->m_traits->parameter_space_in_y_2_object()(overlap_cv, ARR_MAX_END);
  if ((ps_x_r != ARR_INTERIOR) || (ps_y_r != ARR_INTERIOR)) {
    CGAL_assertion(c1->right_event() == c2->right_event());
    right_event = (Event*)(c1->right_event());
  }
  else {
    Point_2 right_end = this->m_traits->construct_max_vertex_2_object()(overlap_cv);
    right_event = this->_push_event(right_end, Event::DEFAULT, ARR_INTERIOR, ARR_INTERIOR).first;
  }

  if (!c1->is_start_point(left_event))
  {
    // here we do not add a curve on the left if there wasn't a curve before
    // it might happen that a curve will be added on the left while
    // it should have been an overlapping curve (that will be detected
    // upon handling of an event with no left curve).
    // See for example data/test_construction/segments/test51.txt
    if (!c2->is_start_point(left_event) || left_event->has_left_curves())
      left_event->add_curve_to_left(c1);
  }
  else
    left_event->remove_curve_from_right(c1);
  if (!c2->is_start_point(left_event))
  {
    // same reason than just above
    if (!c1->is_start_point(left_event) || left_event->has_left_curves())
      left_event->add_curve_to_left(c2);
  }
  else
    left_event->remove_curve_from_right(c2);

  // Allocate the new Subcurve for the overlap
  Subcurve* overlap_sc=NULL;
  if (all_leaves_diff.empty())
  {
    // first check that an equivalent curve is not already in left_event
    for (Subcurve_iterator iter = left_event->right_curves_begin(); iter != left_event->right_curves_end();
     ++iter)
    {
      if ( (*iter)->has_same_leaves(c1, c2) )
      {
        CGAL_SS_PRINT_TEXT("Reuse overlapping curve ");
        CGAL_SS_PRINT_CURVE(*iter);
        CGAL_SS_PRINT_EOL();
        overlap_sc=*iter;
        break;
      }
    }

    if (overlap_sc==NULL)
    {
      CGAL_SS_PRINT_TEXT("Allocate a new subcurve for the overlap (no common subcurves)");
      CGAL_SS_PRINT_EOL();
      // no duplicate only one curve is needed
      overlap_sc = this->m_subCurveAlloc.allocate(1);
#ifdef CGAL_CXX11
      std::allocator_traits<Subcurve_alloc>::construct(this->m_subCurveAlloc,overlap_sc, this->m_masterSubcurve);
#else
      this->m_subCurveAlloc.construct(overlap_sc, this->m_masterSubcurve);
#endif
      overlap_sc->set_hint(this->m_statusLine.end());
      overlap_sc->init(overlap_cv);
      overlap_sc->set_left_event(left_event);
      overlap_sc->set_right_event(right_event);
      m_overlap_subCurves.push_back(overlap_sc);
      // sets the two originating subcurves of overlap_sc
      overlap_sc->set_originating_subcurve1(c1);
      overlap_sc->set_originating_subcurve2(c2);
    }
  }
  else{
    CGAL_SS_PRINT_TEXT("Allocate new subcurves for the overlap (common subcurves)");
    CGAL_SS_PRINT_EOL();

    // create an overlapping curve per subcurve in second_parent that is not in first_parent
    for (typename std::vector<Subcurve*>::const_iterator sc_it=all_leaves_diff.begin();
                                                         sc_it!=all_leaves_diff.end();
                                                         ++sc_it)
    {
      overlap_sc = this->m_subCurveAlloc.allocate(1);
#ifdef CGAL_CXX11
      std::allocator_traits<Subcurve_alloc>::construct(this->m_subCurveAlloc,overlap_sc, this->m_masterSubcurve);
#else
      this->m_subCurveAlloc.construct(overlap_sc, this->m_masterSubcurve);
#endif
      overlap_sc->set_hint(this->m_statusLine.end());
      overlap_sc->init(overlap_cv);
      overlap_sc->set_left_event(left_event);
      overlap_sc->set_right_event(right_event);
      m_overlap_subCurves.push_back(overlap_sc);
      // sets the two originating subcurves of overlap_sc
      overlap_sc->set_originating_subcurve1(first_parent);
      overlap_sc->set_originating_subcurve2(*sc_it);
      first_parent=overlap_sc;
    }
  }
  left_event->set_overlap();

  CGAL_SS_PRINT_CURVE(c1);
  CGAL_SS_PRINT_TEXT(" + ");
  CGAL_SS_PRINT_CURVE(c2);
  CGAL_SS_PRINT_TEXT(" => ");
  CGAL_SS_PRINT_EOL();
  CGAL_SS_PRINT_TEXT("  ");
  CGAL_SS_PRINT_CURVE(overlap_sc);
  CGAL_SS_PRINT_EOL();

  // Remove curves from the left curves of the right end
  // and add them on the right otherwise
  if (c1->is_end_point(right_event))
    right_event->remove_curve_from_left(c1);
  else
    _add_curve_to_right(right_event, c1);

  if (c2->is_end_point(right_event))
    right_event->remove_curve_from_left(c2);
  else
    _add_curve_to_right(right_event, c2);

  // add the overlapping curve of the right of the left end
  _add_curve_to_right(left_event, overlap_sc);
  right_event->add_curve_to_left(overlap_sc);

  this->m_visitor->found_overlap(c1, c2, overlap_sc);

  if (!c1->is_end_point(right_event) && !c2->is_end_point(right_event))
    if (right_event->is_right_curve_bigger(c1, c2, this->m_traits))
      std::swap(c1, c2);

  if (event_on_overlap != NULL && event_on_overlap!=left_event)
  {
    // first check that event_on_overlap is in that overlapping curve as it
    // can be on another overlapping part of the generating curves
    if ( this->m_queueEventLess(event_on_overlap, left_event) == LARGER &&
         this->m_queueEventLess(right_event, event_on_overlap) == LARGER )
    {
      CGAL_SS_PRINT_TEXT("An overlapping curve starts before the event inducing its creation");
      CGAL_SS_PRINT_EOL();
      // an overlapping curve is involved in an intersection
      _add_curve_to_right(event_on_overlap, overlap_sc);
      event_on_overlap->add_curve_to_left(overlap_sc);
    }
  }
  CGAL_SS_PRINT_END_EOL("creating an overlapping curve");
}

//-----------------------------------------------------------------------------
// Fix a subcurve that represents an overlap.
// sc - some originating subcurve of a aubcurve that stores an overlap
// notice thah this function is recursive since an originating subcurve of
// an overlap can be itself a subcurve that stores overlap and so on.
template <typename Vis>
void Surface_sweep_2<Vis>::_fix_finished_overlap_subcurve(Subcurve* sc)
{
  CGAL_SS_PRINT_START("fixing finished overlap subcurve ");
  CGAL_SS_PRINT_CURVE(sc);
  CGAL_SS_PRINT_EOL();

  CGAL_assertion(sc != NULL);

  // split 'sc' if necessary and update to event as weak intersection
  if ((Event*)sc->right_event() != this->m_currentEvent) {
    this->m_traits->split_2_object()(sc->last_curve(),
                                     this->m_currentEvent->point(),
                                     sub_cv1, sub_cv2);
    sc->set_last_curve(sub_cv2);

    this->m_currentEvent->set_weak_intersection();
    this->m_visitor->update_event(this->m_currentEvent, sc);

    CGAL_SS_PRINT_END_EOL("Fixing finished overlap subcurve");
    return;
  }

  if (!sc->originating_subcurve1()) {
    // sc does not store an overlap, we are done
    CGAL_SS_PRINT_END_EOL("fixing finished overlap subcurve");
    return;
  }

  // sc is a subcurve that stores overlap, we have to continue with the
  // recursion and deal with his two originating subcurves recursively.
  Subcurve* orig_sc_1 = sc->originating_subcurve1();
  Subcurve* orig_sc_2 = sc->originating_subcurve2();

  _fix_finished_overlap_subcurve(orig_sc_1);
  _fix_finished_overlap_subcurve(orig_sc_2);

  CGAL_SS_PRINT_END_EOL("fixing finished overlap subcurve");
}

//-----------------------------------------------------------------------------
// add a curve as a right curve or left curve when the event is created
// or updated.
//
template <typename Vis>
void Surface_sweep_2<Vis>::_add_curve(Event* e, Subcurve* sc, Attribute type)
{
  if (sc == NULL) return;

  if (type == Event::LEFT_END) {
    sc->set_left_event(e);
    _add_curve_to_right(e, sc);
    return;
  }

  CGAL_assertion(type == Event::RIGHT_END);
  sc->set_right_event(e);
  e->add_curve_to_left(sc);
}

} // namespace Surface_sweep_2
} // namespace CGAL

#endif
