// Copyright (c) 2007,2009,2010,2011 Tel-Aviv University (Israel).
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
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Ron Wein <wein@post.tau.ac.il>
//                 Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_ARR_SPHERICAL_OVERLAY_HELPER_H
#define CGAL_ARR_SPHERICAL_OVERLAY_HELPER_H

#include <CGAL/license/Arrangement_on_surface_2.h>

/*! \file
 *
 * Definition of the Arr_spherical_overlay_helper class-template.
 */

#include <CGAL/Arr_topology_traits/Arr_spherical_construction_helper.h>

namespace CGAL {

/*! \class Arr_spherical_overlay_helper
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
class Arr_spherical_overlay_helper {
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

  typedef typename Event::Subcurve_iterator             Subcurve_iterator;


  // The input arrangements (the "red" and the "blue" one):
  typedef typename Ar2::Topology_traits                 Topology_traits_red;
  typedef typename Ar2::Face_const_handle               Face_handle_red;

  typedef typename Ab2::Topology_traits                 Topology_traits_blue;
  typedef typename Ab2::Face_const_handle               Face_handle_blue;

  // Define the helper class for the construction visitor.
  typedef Arr_spherical_construction_helper<Gt2, Arrangement_2, Event,
                                            Subcurve>     Construction_helper;

protected:
  // Data members:
  const Topology_traits_red* m_red_top_traits;
  const Topology_traits_blue* m_blue_top_traits;

  //! Red spherical face
  Face_handle_red m_red_nf;

  //! Blue spherical face
  Face_handle_blue m_blue_nf;

public:
  /*! Constructor, given the input red and blue arrangements. */
  Arr_spherical_overlay_helper(const Ar2* red_arr, const Ab2* blue_arr) :
    m_red_top_traits(red_arr->topology_traits()),
    m_blue_top_traits(blue_arr->topology_traits())
  {}

  /// \name Notification functions.
  //@{

  /* A notification issued before the sweep process starts. */
  void before_sweep()
  {
    // Get the spherical faces in both arrangements.
    /* RWRW:
     * m_red_nf = Face_handle_red(m_red_top_traits->spherical_face());
     * m_blue_nf = Face_handle_blue(m_blue_top_traits->spherical_face());
     */
    m_red_nf = Face_handle_red(m_red_top_traits->south_face());
    m_blue_nf = Face_handle_blue(m_blue_top_traits->south_face());
  }

  /*! A notification invoked before the sweep-line starts handling the given
   * event.
   */
  void before_handle_event(Event* event)
  {
    if (event->parameter_space_in_y() != ARR_TOP_BOUNDARY &&
        event->parameter_space_in_x() != ARR_LEFT_BOUNDARY)
      return;

    Arr_curve_end ind = ((event->number_of_left_curves() == 0) &&
                         (event->number_of_right_curves() != 0)) ?
      ARR_MIN_END : ARR_MAX_END;

    Subcurve_iterator it_red, it_blue, it_end;
    if (ind == ARR_MIN_END) {
      it_blue = it_red = event->right_curves_begin();
      it_end = event->right_curves_end();
    }
    else {
      it_blue = it_red = event->left_curves_begin();
      it_end = event->left_curves_end();
    }

    // red arrangement
    while ((it_red != it_end) && ((*it_red)->color() == Gt2::BLUE))
      ++it_red;

    if (it_red != it_end) {
      const Subcurve* sc_red = *it_red;
      if (event->parameter_space_in_y() == ARR_TOP_BOUNDARY) {
        // The curve is incident to the north pole.
        switch (sc_red->color()) {
         case Gt2::RED:
          m_red_nf = (ind == ARR_MIN_END) ?
            sc_red->red_halfedge_handle()->twin()->face() :
            sc_red->red_halfedge_handle()->face();
          break;

         case Gt2::RB_OVERLAP:
          m_red_nf = (ind == ARR_MIN_END) ?
            sc_red->red_halfedge_handle()->twin()->face() :
            sc_red->red_halfedge_handle()->face();
          break;

         case Gt2::BLUE: break;
        }
      }
      else {
        // The curve extends to the right from the curve of discontinuity.
        CGAL_assertion(ind == ARR_MIN_END);
        switch (sc_red->color()) {
        case Gt2::RED:
          m_red_nf = sc_red->red_halfedge_handle()->twin()->face();
          break;
        case Gt2::RB_OVERLAP:
          m_red_nf = sc_red->red_halfedge_handle()->twin()->face();
          break;
        case Gt2::BLUE: break;
        }
      }
    }

    // blue arrangement
    while ((it_blue != it_end) && ((*it_blue)->color() == Gt2::RED))
      ++it_blue;

    if (it_blue != it_end) {
      const Subcurve* sc_blue = *it_blue;
      if (event->parameter_space_in_y() == ARR_TOP_BOUNDARY) {
        // The curve is incident to the north pole.
        switch (sc_blue->color()) {
         case Gt2::BLUE:
          m_blue_nf = (ind == ARR_MIN_END) ?
            sc_blue->blue_halfedge_handle()->twin()->face() :
            sc_blue->blue_halfedge_handle()->face();
          break;

         case Gt2::RB_OVERLAP:
          m_blue_nf = (ind == ARR_MIN_END) ?
            sc_blue->blue_halfedge_handle()->twin()->face() :
            sc_blue->blue_halfedge_handle()->face();
          break;

         case Gt2::RED: break;
        }
      }
      else {
        // The curve extends to the right from the curve of discontinuity.
        CGAL_assertion(ind == ARR_MIN_END);
        switch (sc_blue->color()) {
         case Gt2::BLUE:
          m_blue_nf = sc_blue->blue_halfedge_handle()->twin()->face();
          break;

         case Gt2::RB_OVERLAP:
          m_blue_nf = sc_blue->blue_halfedge_handle()->twin()->face();
          break;

         case Gt2::RED: break;
        }
      }
    }
  }
  //@}

  /*! Get the current red top face. */
  Face_handle_red red_top_face() const { return m_red_nf; }

  /*! Get the current blue top face. */
  Face_handle_blue blue_top_face() const { return m_blue_nf; }

  /*! Obtain the red topology traits. */
  const Topology_traits_red* red_topology_traits() const
  { return m_red_top_traits; }

  /*! Obtain the blue topology traits. */
  const Topology_traits_blue* blue_topology_traits() const
  { return m_blue_top_traits; }
};

} //namespace CGAL

#endif
