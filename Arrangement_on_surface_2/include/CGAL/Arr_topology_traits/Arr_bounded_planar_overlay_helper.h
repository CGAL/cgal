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
//                 Ron Wein <wein@post.tau.ac.il>
//                 Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_ARR_BOUNDED_PLANAR_OVERLAY_HELPER_H
#define CGAL_ARR_BOUNDED_PLANAR_OVERLAY_HELPER_H

#include <CGAL/license/Arrangement_on_surface_2.h>

/*! \file
 *
 * Definition of the Arr_bounded_planar_overlay_helper class-template.
 */

#include <CGAL/Arr_topology_traits/Arr_bounded_planar_construction_helper.h>

namespace CGAL {

/*! \class Arr_bounded_planar_overlay_helper
 *
 * A helper class for the overlay sweep-line visitor, suitable for the overlay
 * of Arrangement_on_surface_2 objects instantiated with a topology-traits
 * class for bounded curves in the plane.
 */
template <typename GeometryTraits_2,
          typename ArrangementRed_2,
          typename ArrangementBlue_2,
          typename Arrangement_,
          typename Event_,
          typename Subcurve_>
class Arr_bounded_planar_overlay_helper {
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
  typedef typename Ar2::Face_const_handle               Face_handle_red;

  typedef typename Ab2::Face_const_handle               Face_handle_blue;

  // Define the helper class for the construction visitor.
  typedef Arr_bounded_planar_construction_helper<Gt2, Arrangement_2, Event,
                                                 Subcurve>
                                                        Construction_helper;

protected:
  // Data members:
  const typename Ar2::Topology_traits* m_red_top_traits;
  const typename Ab2::Topology_traits* m_blue_top_traits;

  Face_handle_red m_red_ubf;            // Red unbounded face.
  Face_handle_blue m_blue_ubf;          // Blue unbounded face.

public:
  /*! Constructor, given the input red and blue arrangements. */
  Arr_bounded_planar_overlay_helper(const Ar2* red_arr, const Ab2* blue_arr) :
    m_red_top_traits (red_arr->topology_traits()),
    m_blue_top_traits (blue_arr->topology_traits())
  {}

  /// \name Notification functions.
  //@{

  /* A notification issued before the sweep process starts. */
  void before_sweep()
  {
    // Get the unbounded faces in both arrangements.
    m_red_ubf = Face_handle_red (m_red_top_traits->unbounded_face());
    m_blue_ubf = Face_handle_blue (m_blue_top_traits->unbounded_face());
  }

  /*! A notification invoked before the sweep-line starts handling the given
   * event.
   */
  void before_handle_event(Event* /* e */) {}
  //@}

  /*! Obtain the current red top face. */
  Face_handle_red red_top_face() const { return m_red_ubf; }

  /*! Obtain the current blue top face. */
  Face_handle_blue blue_top_face() const { return m_blue_ubf; }
};

} // namespace CGAL

#endif
