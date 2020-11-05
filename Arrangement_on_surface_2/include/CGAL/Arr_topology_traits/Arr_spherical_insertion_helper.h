// Copyright (c) 2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Ron Wein <wein@post.tau.ac.il>
//                 Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_ARR_SPHERICAL_INSERTION_HELPER_H
#define CGAL_ARR_SPHERICAL_INSERTION_HELPER_H

#include <CGAL/license/Arrangement_on_surface_2.h>

/*! \file
 *
 * Definition of the Arr_spherical_insertion_helper class-template.
 */

#include <CGAL/Arr_topology_traits/Arr_spherical_construction_helper.h>

namespace CGAL {

/*! \class Arr_spherical_insertion_helper
 *
 * A helper class for the insertion sweep-line visitors, suitable
 * for an Arrangement_on_surface_2 instantiated with a topology-traits class
 * for bounded curves in the plane.
 */
template <typename GeometryTraits_2, typename Arrangement_, typename Event_,
          typename Subcurve_>
class Arr_spherical_insertion_helper :
  public Arr_spherical_construction_helper<GeometryTraits_2, Arrangement_,
                                           Event_, Subcurve_>
{
public:
  typedef GeometryTraits_2                              Geometry_traits_2;
  typedef Arrangement_                                  Arrangement_2;
  typedef Event_                                        Event;
  typedef Subcurve_                                     Subcurve;

private:
  typedef Geometry_traits_2                             Gt2;
  typedef Arr_spherical_insertion_helper<Gt2, Arrangement_2, Event, Subcurve>
                                                        Self;
  typedef Arr_spherical_construction_helper<Gt2, Arrangement_2, Event, Subcurve>
                                                        Base;

public:
  typedef typename Gt2::X_monotone_curve_2              X_monotone_curve_2;
  typedef typename Gt2::Point_2                         Point_2;

  typedef typename Arrangement_2::Vertex_handle         Vertex_handle;
  typedef typename Arrangement_2::Face_handle           Face_handle;
  typedef typename Arrangement_2::Halfedge_handle       Halfedge_handle;

protected:
  typedef typename Arrangement_2::Topology_traits       Topology_traits;
  typedef typename Topology_traits::Vertex              DVertex;

  //! The halfedge that points at the spherical face.
  Halfedge_handle m_spherical_halfedge;

public:
  /*! Constructor */
  Arr_spherical_insertion_helper(Arrangement_2 *arr) : Base(arr) {}

  /*! Destructor. */
  virtual ~Arr_spherical_insertion_helper() {}

  /// \name Notification functions.
  //@{

  /* A notification issued before the sweep process starts. */
  virtual void before_sweep()
  {
    // Get the unbounded face.
    this->m_spherical_face = Face_handle(this->m_top_traits->south_face());
  }

  /*! A notification invoked before the sweep-line starts handling a given
   * event.
   */
  virtual void before_handle_event(Event* event);

  //@}

  /*! A notification invoked when a new subcurve is created. */
  virtual void add_subcurve(Halfedge_handle he, Subcurve* sc);

  /*! Get the current top face. */
  virtual Face_handle top_face() const;
};

//-----------------------------------------------------------------------------
// Memeber-function definitions:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// A notification invoked before the sweep-line starts handling the given
// event.
//
template <typename Tr, typename Arr, typename Evnt, typename Sbcv>
void Arr_spherical_insertion_helper<Tr, Arr, Evnt, Sbcv>::
before_handle_event(Event* event)
{
  // Ignore events that do not have boundary conditions.
  const Arr_parameter_space ps_x = event->parameter_space_in_x();
  const Arr_parameter_space ps_y = event->parameter_space_in_y();
  if ((ps_x == ARR_INTERIOR) && (ps_y == ARR_INTERIOR)) return;

  if (ps_y == ARR_BOTTOM_BOUNDARY) {
    // Process bootom contraction boundary:
    // The event has only one right curve, as there is exactly one curve
    // incident to an event with boundary conditions.
    CGAL_assertion((event->number_of_left_curves() == 0) &&
                   (event->number_of_right_curves() == 1));
    const X_monotone_curve_2& xc =
        (*(event->right_curves_begin()))->last_curve();
    if (xc.halfedge_handle() != Halfedge_handle()) return;

    // If a vertex on the south pole does not exists, create one.
    DVertex* dv = this->m_top_traits->south_pole();
    Vertex_handle v = (dv) ? Vertex_handle(dv) :
      this->m_arr_access.create_boundary_vertex(xc, ARR_MIN_END, ps_x, ps_y);
    event->set_vertex_handle(v);
    return;
  }

  if (ps_y == ARR_TOP_BOUNDARY) {
    // Process top contraction boundary:
    // The event has only one left curve, as there is exactly one curve
    // incident to an event with boundary conditions.
    CGAL_assertion((event->number_of_left_curves() == 1) &&
                   (event->number_of_right_curves() == 0));
    const X_monotone_curve_2& xc =
      (*(event->left_curves_begin()))->last_curve();

    if (xc.halfedge_handle() != Halfedge_handle()) {
      // Update the current top face.
      this->m_spherical_face = xc.halfedge_handle()->face();
      return;
    }

    // If a vertex on the north pole does not exist, create one.
    DVertex* dv = this->m_top_traits->north_pole();
    Vertex_handle v = (dv) ? Vertex_handle(dv) :
      this->m_arr_access.create_boundary_vertex(xc, ARR_MAX_END, ps_x, ps_y);
    event->set_vertex_handle(v);
    return;
  }

  if (ps_x == ARR_LEFT_BOUNDARY) {
    // Process left discontinuity boundary:
    // The event has only right curves, as there is exactly one curve
    // incident to an event with boundary conditions.
    CGAL_assertion((event->number_of_left_curves() == 0) &&
                   (event->number_of_right_curves() >= 1));
    const X_monotone_curve_2& xc =
      (*(event->right_curves_begin()))->last_curve();

    if (xc.halfedge_handle() != Halfedge_handle()) {
      // Update the current top face.
      this->m_spherical_face = xc.halfedge_handle()->twin()->face();
      return;
    }

    // If a vertex on the line of discontinuity does not exists. create one.
    DVertex* dv = this->m_top_traits->discontinuity_vertex(xc, ARR_MIN_END);
    Vertex_handle v = (dv) ? Vertex_handle(dv) :
      this->m_arr_access.create_boundary_vertex(xc, ARR_MIN_END, ps_x, ps_y);
    event->set_vertex_handle(v);
    return;
  }

  if (ps_x == ARR_RIGHT_BOUNDARY) {
    // Process right discontinuity boundary:
    // The event has only left curves, as there is exactly one curve
    // incident to an event with boundary conditions.
    CGAL_assertion((event->number_of_left_curves() >= 1) &&
                   (event->number_of_right_curves() == 0));
    const X_monotone_curve_2& xc =
      (*(event->left_curves_begin()))->last_curve();
    if (xc.halfedge_handle() != Halfedge_handle()) return;

    // If a vertex on the line of discontinuity does not exists. create one.
    DVertex* dv = this->m_top_traits->discontinuity_vertex(xc, ARR_MAX_END);
    Vertex_handle v = (dv) ? Vertex_handle(dv) :
      this->m_arr_access.create_boundary_vertex(xc, ARR_MAX_END, ps_x, ps_y);
    event->set_vertex_handle(v);
    return;
  }
}

/*! A notification invoked when a new subcurve is created. */
template <typename Tr, typename Arr, typename Evnt, typename Sbcv>
void Arr_spherical_insertion_helper<Tr, Arr, Evnt, Sbcv>::
add_subcurve(Halfedge_handle he, Subcurve* /* sc */)
{
  if (he->source()->parameter_space_in_y() == ARR_TOP_BOUNDARY) {
    m_spherical_halfedge = he;
    return;
  }
  if (he->target()->parameter_space_in_y() == ARR_TOP_BOUNDARY) {
    m_spherical_halfedge = he->twin();
    return;
  }
}

/*! Get the current top face. */
template <typename Tr, typename Arr, typename Evnt, typename Sbcv>
typename Arr_spherical_insertion_helper<Tr, Arr, Evnt, Sbcv>::Face_handle
Arr_spherical_insertion_helper<Tr, Arr, Evnt, Sbcv>::top_face() const
{
  const Halfedge_handle invalid_he;
  if (m_spherical_halfedge != invalid_he) return m_spherical_halfedge->face();
  return this->m_spherical_face;
}

} // namespace CGAL

#endif
