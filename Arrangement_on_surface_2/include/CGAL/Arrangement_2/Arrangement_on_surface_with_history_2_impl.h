// Copyright (c) 2005,2007,2008,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s): Ron Wein         <wein@post.tau.ac.il>
//            Efi Fogel        <efif@post.tau.ac.il>
//            Baruch Zukerman  <baruchzu@post.tau.ac.il>

#ifndef CGAL_ARRANGEMENT_ON_SURFACE_WITH_HISTORY_2_IMPL_H
#define CGAL_ARRANGEMENT_ON_SURFACE_WITH_HISTORY_2_IMPL_H

#include <CGAL/license/Arrangement_on_surface_2.h>

/*! \file
 * Member-function definitions for the `Arrangement_on_surface_with_history_2`
 * class.
 */

namespace CGAL {

//-----------------------------------------------------------------------------
// Default constructor.
//
template <typename GeomTraits, typename TopolTraits>
Arrangement_on_surface_with_history_2<GeomTraits, TopolTraits>::Arrangement_on_surface_with_history_2()
{ m_observer.attach(*this); }

//-----------------------------------------------------------------------------
// Copy constructor.
//
template <typename GeomTraits, typename TopolTraits>
Arrangement_on_surface_with_history_2<GeomTraits, TopolTraits>::Arrangement_on_surface_with_history_2(const Self& arr) :
  Base_arr_2(arr.Base_arr_2::shared_geometry_traits()) {
  assign(arr);
  m_observer.attach(*this);
}

//-----------------------------------------------------------------------------
// Constructor given a shared traits object.
//
template <typename GeomTraits, typename TopolTraits>
Arrangement_on_surface_with_history_2<GeomTraits, TopolTraits>::
Arrangement_on_surface_with_history_2(Shared_geometry_traits tr) :
  Base_arr_2(_data_traits(tr))
{ m_observer.attach(*this); }

//-----------------------------------------------------------------------------
// Constructor given a traits object owned by the caller.
//
template <typename GeomTraits, typename TopolTraits>
Arrangement_on_surface_with_history_2<GeomTraits, TopolTraits>::
Arrangement_on_surface_with_history_2(const Geometry_traits_2* tr) :
  Base_arr_2(static_cast<const Data_traits_2*>(tr))
{ m_observer.attach(*this); }

//-----------------------------------------------------------------------------
// Assignment operator.
//
template <typename GeomTraits, typename TopolTraits>
Arrangement_on_surface_with_history_2<GeomTraits, TopolTraits>&
Arrangement_on_surface_with_history_2<GeomTraits, TopolTraits>::operator=(const Self& arr) {
  if (this != &arr) assign(arr);
  return *this;
}

//-----------------------------------------------------------------------------
// Assign an arrangement with history.
//
template <typename GeomTraits, typename TopolTraits>
void Arrangement_on_surface_with_history_2<GeomTraits, TopolTraits>::assign(const Self& arr) {
  // Clear the current contents of the arrangement.
  clear();

  // Assign the base arrangement.
  Base_arr_2::assign(arr);

  // Duplicate the curves of arr, and redirect the curve pointers stored with the edges (which were copied from arr)
  // to the duplicates.
  Curve_map cv_map;
  _duplicate_curves(arr.curves_begin(), arr.curves_end(), cv_map);
  _relink_edges(cv_map);
}

//-----------------------------------------------------------------------------
// Destructor.
//
template <typename GeomTraits, typename TopolTraits>
Arrangement_on_surface_with_history_2<GeomTraits, TopolTraits>::~Arrangement_on_surface_with_history_2() { clear(); }

//-----------------------------------------------------------------------------
// Clear the arrangement.
//
template <typename GeomTraits, typename TopolTraits>
void Arrangement_on_surface_with_history_2<GeomTraits, TopolTraits>::clear() {
  // Free all stored curves. Note that the iterator is advanced before the curve it points to is erased.
  for (auto cit = m_curves.begin(); cit != m_curves.end();) _delete_curve_halfedges(&*cit++);

  // Clear the base arrangement.
  Base_arr_2::clear();
}

//-----------------------------------------------------------------------------
// Split a given edge into two at the given split point.
//
template <typename GeomTraits, typename TopolTraits>
typename Arrangement_on_surface_with_history_2<GeomTraits, TopolTraits>::Halfedge_handle
Arrangement_on_surface_with_history_2<GeomTraits, TopolTraits>::split_edge(Halfedge_handle e, const Point_2& p) {
  // Split the curve associated with the halfedge e at the given point p.
  Data_x_curve_2 cv1;
  Data_x_curve_2 cv2;

  this->m_geom_traits->split_2_object()(e->curve(), p, cv1, cv2);

  // cv1 always lies to the left of cv2. If e is directed from left to right,
  // we should split and return the halfedge associated with cv1, and
  // otherwise we should return the halfedge associated with cv2 after the
  // split.
  return (e->direction() == ARR_LEFT_TO_RIGHT) ?
    Base_arr_2::split_edge(e, cv1, cv2) : Base_arr_2::split_edge(e, cv2, cv1);
}

//-----------------------------------------------------------------------------
// Merge two edges to form a single edge.
//
template <typename GeomTraits, typename TopolTraits>
typename Arrangement_on_surface_with_history_2<GeomTraits, TopolTraits>::Halfedge_handle
Arrangement_on_surface_with_history_2<GeomTraits, TopolTraits>::merge_edge(Halfedge_handle e1, Halfedge_handle e2) {
  CGAL_precondition_msg(are_mergeable(e1, e2), "Edges are not mergeable.");

  // Merge the two curves.
  Data_x_curve_2 cv;
  this->m_geom_traits->merge_2_object()(e1->curve(), e2->curve(), cv);
  return Base_arr_2::merge_edge(e1, e2, cv);
}

//-----------------------------------------------------------------------------
// Check if two edges can be merged to a single edge.
//
template <typename GeomTraits, typename TopolTraits>
bool Arrangement_on_surface_with_history_2<GeomTraits, TopolTraits>::
are_mergeable(Halfedge_const_handle e1, Halfedge_const_handle e2) const {
  // Both halfedges must be non-fictitious.
  if (e1->is_fictitious() || e2->is_fictitious()) return false;

  // In order to be mergeable, the two halfedges must share a common
  // end-vertex. We assign vh to be this vertex.
  Vertex_const_handle vh;

  if (e1->target() == e2->source() || e1->target() == e2->target()) vh = e1->target();
  else if (e1->source() == e2->source() || e1->source() == e2->target()) vh = e1->source();
  else return false;  // No common end-vertex: the edges are not mergeable.

  // If there are other edges incident to vh, it is impossible to remove it
  // and merge the two edges.
  if (vh->degree() != 2) return false;

  // Check whether the curves associated with the two edges are mergeable.
  return this->m_geom_traits->are_mergeable_2_object()(e1->curve(), e2->curve());
}

} // namespace CGAL

#endif
