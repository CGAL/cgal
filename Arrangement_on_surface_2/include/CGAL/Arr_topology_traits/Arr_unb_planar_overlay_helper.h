// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Ron Wein        <wein@post.tau.ac.il>
//                 Efi Fogel       <efif@post.tau.ac.il>

#ifndef CGAL_ARR_UNB_PLANAR_OVERLAY_HELPER_H
#define CGAL_ARR_UNB_PLANAR_OVERLAY_HELPER_H

#include <CGAL/license/Arrangement_on_surface_2.h>

/*! \file
 *
 * Definition of the Arr_unb_planar_overlay_helper class-template.
 */

#include <CGAL/Arr_topology_traits/Arr_unb_planar_construction_helper.h>

namespace CGAL {

/*! \class Arr_unb_planar_overlay_helper
 *
 * A helper class for the overlay sweep-line visitor, suitable for the overlay
 * of Arrangement_on_surface_2 objects instantiated with a topology-traits
 * class for unbounded curves in the plane.
 */
template <typename GeometryTraits_2,
          typename ArrangementRed_2,
          typename ArrangementBlue_2,
          typename Arrangement_,
          typename Event_,
          typename Subcurve_>
class Arr_unb_planar_overlay_helper {
public:
  typedef GeometryTraits_2                              Geometry_traits_2;
  typedef ArrangementRed_2                              Arrangement_red_2;
  typedef ArrangementBlue_2                             Arrangement_blue_2;
  typedef Arrangement_                                  Arrangement_2;
  typedef Event_                                        Event;
  typedef Subcurve_                                     Subcurve;

private:
  typedef Geometry_traits_2                             Gt2;
  typedef Arrangement_red_2                             Ar2;
  typedef Arrangement_blue_2                            Ab2;

public:
  typedef typename Gt2::X_monotone_curve_2              X_monotone_curve_2;
  typedef typename Gt2::Point_2                         Point_2;

  // The input arrangements (the "red" and the "blue" one):
  typedef typename Ar2::Halfedge_const_handle           Halfedge_handle_red;
  typedef typename Ar2::Face_const_handle               Face_handle_red;
  typedef typename Ar2::Vertex_const_handle             Vertex_handle_red;

  typedef typename Ab2::Halfedge_const_handle           Halfedge_handle_blue;
  typedef typename Ab2::Face_const_handle               Face_handle_blue;
  typedef typename Ab2::Vertex_const_handle             Vertex_handle_blue;

  // Define the helper class for the construction visitor.
  typedef Arr_unb_planar_construction_helper<Gt2, Arrangement_2, Event,
                                             Subcurve>  Construction_helper;

protected:
  // Data members:
  const typename Ar2::Topology_traits* m_red_top_traits;
  const typename Ab2::Topology_traits* m_blue_top_traits;

  Halfedge_handle_red m_red_th;         // Red top fictitious halfedge.
  Halfedge_handle_blue m_blue_th;       // Blue top fictitious halfedge.

  Vertex_handle_red v_red_tl;           // Red top-left fictitious vertex.
  Vertex_handle_blue v_blue_tl;         // Blue top-left fictitious vertex.

public:
  /*! Constructor, given the input red and blue arrangements. */
  Arr_unb_planar_overlay_helper(const Ar2* red_arr, const Ab2* blue_arr) :
    m_red_top_traits(red_arr->topology_traits()),
    m_blue_top_traits(blue_arr->topology_traits())
  {}

  /// \name Notification functions.
  //@{

  /* A notification issued before the sweep process starts. */
  void before_sweep();

  /*! A notification invoked before the sweep-line starts handling the given
   * event.
   */
  void before_handle_event(Event* e);
  //@}

  /*! Get the current red top face. */
  Face_handle_red red_top_face() const { return (m_red_th->face()); }

  /*! Get the current blue top face. */
  Face_handle_blue blue_top_face() const { return (m_blue_th->face()); }
};

//-----------------------------------------------------------------------------
// Member-function definitions:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// A notification issued before the sweep process starts.
//
template <typename Tr, typename ArrR, typename ArrB, typename Arr,
          typename Evnt, typename Sbcv>
void Arr_unb_planar_overlay_helper<Tr, ArrR, ArrB, Arr, Evnt, Sbcv>::
before_sweep()
{
  // Get the top-left fictitious vertices in both arrangements.
  v_red_tl = Vertex_handle_red(m_red_top_traits->top_left_vertex());
  v_blue_tl = Vertex_handle_blue(m_blue_top_traits->top_left_vertex());

  // Get the fictitious halfedges incident to the bottom-left fictitious
  // vertices in both red and blue arrangements. If there are no vertices
  // at x = -oo, we take the halfedge incident to the top-left vertex that
  // lies on the top edge of the fictitious face.
  Vertex_handle_red v_red_bl =
    Vertex_handle_red(m_red_top_traits->bottom_left_vertex());

  m_red_th = v_red_bl->incident_halfedges();

  if (m_red_th->source()->parameter_space_in_x() != ARR_LEFT_BOUNDARY)
    m_red_th = m_red_th->next()->twin();

  if (m_red_th->source() == v_red_tl) m_red_th = m_red_th->prev();

  Vertex_handle_blue v_blue_bl =
    Vertex_handle_blue(m_blue_top_traits->bottom_left_vertex());

  m_blue_th = v_blue_bl->incident_halfedges();
  if (m_blue_th->source()->parameter_space_in_x() != ARR_LEFT_BOUNDARY)
    m_blue_th = m_blue_th->next()->twin();

  if (m_blue_th->source() == v_blue_tl) m_blue_th = m_blue_th->prev();
}

//-----------------------------------------------------------------------------
// A notification invoked before the sweep-line starts handling the given
// event.
//
template <typename Tr, typename ArrR, typename ArrB, typename Arr,
          typename Evnt, typename Sbcv>
void Arr_unb_planar_overlay_helper<Tr, ArrR, ArrB, Arr, Evnt, Sbcv>::
before_handle_event(Event* e)
{
  // Nothing to do in case the event represents a valid point.
  if (e->is_closed()) return;

  // In case the event occurs on the left edge of the fictitious face (x = -oo)
  // or on its top edge (finite x and y = +oo), update the fictitious top
  // halfedges.
  if (e->parameter_space_in_x() == ARR_LEFT_BOUNDARY ||
      (e->parameter_space_in_x() == ARR_INTERIOR &&
       e->parameter_space_in_y() == ARR_TOP_BOUNDARY))
  {
    Arr_curve_end ce;
    switch (e->boundary_touching_curve(ce).color()) {
    case (Gt2::RED) :
      // Update the red top fictitious halfedge.
      m_red_th = m_red_th->twin()->next()->twin();
      if (m_red_th->source() == v_red_tl) m_red_th = m_red_th->prev();
      break;

    case (Gt2::BLUE) :
      // Update the blue top fictitious halfedge.
      m_blue_th = m_blue_th->twin()->next()->twin();
      if (m_blue_th->source() == v_blue_tl) m_blue_th = m_blue_th->prev();
      break;

    case Gt2::RB_OVERLAP :
      // Update both red and blue top fictitious halfedges.
      m_red_th = m_red_th->twin()->next()->twin();
      if (m_red_th->source() == v_red_tl) m_red_th = m_red_th->prev();

      m_blue_th = m_blue_th->twin()->next()->twin();
      if (m_blue_th->source() == v_blue_tl) m_blue_th = m_blue_th->prev();
      break;
    }
  }
}

} // namespace CGAL

#endif
