// Copyright (c) 2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Ron Wein <wein@post.tau.ac.il>
//            Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_ARR_SPHERICAL_VERT_DECOMP_HELPER_H
#define CGAL_ARR_SPHERICAL_VERT_DECOMP_HELPER_H

#include <CGAL/license/Arrangement_on_surface_2.h>

/*! \file
 *
 * Definition of the Arr_spherical_vert_decomp_helper class-template.
 */

namespace CGAL {

/*! \class Arr_spherical_vert_decomp_helper
 *
 * A helper class for the vertical decomposition sweep-line visitor, suitable
 * for an Arrangement_on_surface_2 instantiated with a topology-traits class
 * for bounded curves in the plane.
 */
template <typename GeometryTraits_2, typename Arrangement_, typename Event_,
          typename Subcurve_>
class Arr_spherical_vert_decomp_helper {
public:
  typedef GeometryTraits_2                              Geometry_traits_2;
  typedef Arrangement_                                  Arrangement_2;
  typedef Event_                                        Event;
  typedef Subcurve_                                     Subcurve;
  typedef typename Subcurve::Allocator                  Allocator;

private:
  typedef Geometry_traits_2                             Gt2;

public:
  typedef typename Gt2::X_monotone_curve_2              X_monotone_curve_2;

  typedef typename Arrangement_2::Vertex_const_handle   Vertex_const_handle;
  typedef typename Arrangement_2::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Arrangement_2::Face_const_handle     Face_const_handle;

  typedef boost::variant<Vertex_const_handle, Halfedge_const_handle,
                         Face_const_handle>             Cell_type;
  typedef boost::optional<Cell_type>                    Vert_type;

protected:
  typedef typename Arrangement_2::Topology_traits       Topology_traits;

  const Topology_traits* m_top_traits;        // The topology traits.
  Vertex_const_handle m_north_pole;           // The north pole.
  bool m_valid_north_pole;                    // Is this a valid vertex.
  Face_const_handle m_north_face;             // Current north face.
  Vertex_const_handle m_south_pole;           // The south pole.
  bool m_valid_south_pole;                    // Is this a valid vertex.
  Face_const_handle m_south_face;             // Current south face.

public:
  /*! Constructor.
   * \param arr The arrangement.
   */
  Arr_spherical_vert_decomp_helper(const Arrangement_2* arr) :
    m_top_traits(arr->topology_traits())
  {}

  /// \name Notification functions.
  //@{

  /*! A notification issued before the sweep process starts. */
  void before_sweep();

  /*!
   * A notification invoked after the sweep-line finishes handling the given
   * event.
   */
  void after_handle_event(Event * event);
  //@}

  /*! Get the current top object. */
  Vert_type top_object () const
  {
    return (m_valid_north_pole) ?
      Vert_type(m_north_pole) : Vert_type(m_north_face);
  }

  /*! Get the current bottom object. */
  Vert_type bottom_object () const
  {
    return (m_valid_south_pole) ?
      Vert_type(m_south_pole) : Vert_type(m_south_face);
  }
};

//-----------------------------------------------------------------------------
// Memeber-function definitions:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// A notification issued before the sweep process starts.
//
  template <typename Tr, typename Arr, typename Evnt, typename Sbcv>
  void Arr_spherical_vert_decomp_helper<Tr, Arr, Evnt, Sbcv>::before_sweep()
{
  // Get the north pole and the face that intially contains it.
  m_valid_north_pole = (m_top_traits->north_pole() != nullptr);
  if (m_valid_north_pole)
    m_north_pole = Vertex_const_handle (m_top_traits->north_pole());

  m_north_face = Face_const_handle (m_top_traits->spherical_face());

  // Get the south pole and the face that intially contains it.
  m_valid_south_pole = (m_top_traits->south_pole() != nullptr);
  if (m_valid_south_pole)
    m_south_pole = Vertex_const_handle (m_top_traits->south_pole());

  m_south_face = Face_const_handle (m_top_traits->south_face());
}

//-----------------------------------------------------------------------------
// A notification invoked after the sweep-line finishes handling the given
// event.
///
template <typename Tr, typename Arr, typename Evnt, typename Sbcv>
void Arr_spherical_vert_decomp_helper<Tr, Arr, Evnt, Sbcv>::
after_handle_event(Event *event)
{
  // Ignore events that are not incident to the poles.
  if (event->parameter_space_in_y() == ARR_INTERIOR) return;

  // The is exactly one curve incident to an event with boundary conditions.
  // Obtain this curve and check whether it already exists in the arrangement.
  CGAL_assertion(((event->number_of_left_curves() == 0) &&
                  (event->number_of_right_curves() == 1)) ||
                 ((event->number_of_left_curves() == 1) &&
                  (event->number_of_right_curves() == 0)));

  const Arr_curve_end ind =
    (event->number_of_left_curves() == 0 &&
     event->number_of_right_curves() == 1) ? ARR_MIN_END : ARR_MAX_END;
  const X_monotone_curve_2& xc = (ind == ARR_MIN_END) ?
    (*(event->right_curves_begin()))->last_curve() :
    (*(event->left_curves_begin()))->last_curve();

  if (event->parameter_space_in_y() == ARR_TOP_BOUNDARY) {
    // The event is incident to the north pole: update the north face.
    if (ind == ARR_MIN_END)
      m_north_face = xc.halfedge_handle()->twin()->face();
    else
      m_north_face = xc.halfedge_handle()->face();
  }
  else if (event->parameter_space_in_y() == ARR_BOTTOM_BOUNDARY) {
    // The event is incident to the south pole: update the south face.
    if (ind == ARR_MIN_END) m_south_face = xc.halfedge_handle()->face();
    else m_south_face = xc.halfedge_handle()->twin()->face();
  }
}

} // namespace CGAL

#endif
