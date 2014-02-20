// Copyright (c) 2007,2009,2010,2011,2013,2014 Max-Planck-Institute Saarbruecken (Germany), Tel-Aviv University (Israel).
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
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Ron Wein <wein@post.tau.ac.il>
//                 Efi Fogel <efif@post.tau.ac.il>
//                 Eric Berberich <eric.berberich@cgal.org>

#ifndef CGAL_ARR_TOROIDAL_OVERLAY_HELPER_H
#define CGAL_ARR_TOROIDAL_OVERLAY_HELPER_H

/*! \file
 * Definition of the Arr_toroidal_overlay_helper class-template.
 */

#include <CGAL/Arr_topology_traits/Arr_toroidal_construction_helper.h>

namespace CGAL {

/*! \class Arr_toroidal_overlay_helper
 * A helper class for the overlay sweep-line visitor, suitable for the overlay
 * of Arrangement_on_surface_2 objects instantiated with a topology-traits
 * class for bounded curves in the plane.
 */
template <class Traits_,
          class ArrangementRed_,
          class ArrangementBlue_,
          class Arrangement_,
          class Event_,
          class Subcurve_>
class Arr_toroidal_overlay_helper
{
public:
  typedef Traits_                                         Traits_2;
  typedef Arrangement_                                    Arrangement_2;
  typedef Event_                                          Event;
  typedef Subcurve_                                       Subcurve;

  typedef typename Traits_2::X_monotone_curve_2           X_monotone_curve_2;
  typedef typename Traits_2::Point_2                      Point_2;

  typedef typename Event::Subcurve_iterator               Subcurve_iterator;


  // The input arrangements (the "red" and the "blue" one):
  typedef ArrangementRed_                                 Arrangement_red_2;
  typedef typename Arrangement_red_2::Topology_traits     Topology_traits_red;
  typedef typename Arrangement_red_2::Face_const_handle   Face_handle_red;

  typedef ArrangementBlue_                                Arrangement_blue_2;
  typedef typename Arrangement_blue_2::Topology_traits    Topology_traits_blue;
  typedef typename Arrangement_blue_2::Face_const_handle  Face_handle_blue;

  // Define the helper class for the construction visitor.
  typedef Arr_toroidal_construction_helper<Traits_2, Arrangement_2, Event,
                                            Subcurve>     Construction_helper;

protected:
  // Data members:
  const Topology_traits_red* m_red_top_traits;
  const Topology_traits_blue* m_blue_top_traits;

  //! Red toroidal face
  Face_handle_red m_red_tf;

  //! Blue toroidal face
  Face_handle_blue m_blue_tf;

public:
  /*! Constructor, given the input red and blue arrangements. */
  Arr_toroidal_overlay_helper(const Arrangement_red_2* red_arr,
                               const Arrangement_blue_2* blue_arr) :
    m_red_top_traits(red_arr->topology_traits()),
    m_blue_top_traits(blue_arr->topology_traits())
  {}

  /// \name Notification functions.
  //@{

  /* A notification issued before the sweep process starts. */
  void before_sweep()
  {
  }

  /*! A notification invoked before the sweep-line starts handling the given
   * event.
   */
  void before_handle_event(Event* event)
  {
  }
  //@}

  /*! Get the current red top face. */
  Face_handle_red red_top_face() const { return m_red_tf; }

  /*! Get the current blue top face. */
  Face_handle_blue blue_top_face() const { return m_blue_tf; }

  /*! Obtain the red topology traits. */
  const Topology_traits_red* red_topology_traits() const
  { return m_red_top_traits; }

  /*! Obtain the blue topology traits. */
  const Topology_traits_blue* blue_topology_traits() const
  { return m_blue_top_traits; }

}; // Arr_toroidal_overlay_helper

} //namespace CGAL

#endif // CGAL_ARR_TOROIDAL_OVERLAY_HELPER_H
