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
// Author(s): Baruch Zukerman <baruchzu@post.tau.ac.il>
//            Ron Wein <wein@post.tau.ac.il>
//            Efi Fogel <efif@post.tau.ac.il>

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
  typedef typename Event::Subcurve_reverse_iterator     Subcurve_reverse_iterator;


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

    // 1. The left halfedges and the right halfedges are always directed right
    //    to left.
    // 2. A left curve of an event, if exists, must be vertical, that is, either
    //    a. the event coincides with the top boundary (north pole), or
    //    b. the curve lies on the identification curve
    // 3. The right curves are ordered bottom to top
    // std::cout << "before_handle_event: " << event->point() << std::endl;
    // std::cout << "# left: " << event->number_of_left_curves() << std::endl;
    // for (auto it = event->left_curves_begin();
    //      it != event->left_curves_end(); ++it) {
    //   if ((*it)->color() == Gt2::RED) {
    //     const Subcurve* sc_red = *it;
    //     std::cout << "  Red: "
    //               << sc_red->red_halfedge_handle()->source()->point() << " => "
    //               << sc_red->red_halfedge_handle()->target()->point() << ", "
    //               << sc_red->red_halfedge_handle()->direction()
    //               << std::endl;
    //   }
    //   else {
    //     const Subcurve* sc_blue = *it;
    //     std::cout << "  Blue: "
    //               << sc_blue->blue_halfedge_handle()->source()->point() << " => "
    //               << sc_blue->blue_halfedge_handle()->target()->point() << ", "
    //               << sc_blue->blue_halfedge_handle()->direction()
    //               << std::endl;
    //   }
    // }
    // std::cout << "# right: " << event->number_of_right_curves() << std::endl;
    // for (auto it = event->right_curves_begin();
    //      it != event->right_curves_end(); ++it) {
    //   if ((*it)->color() == Gt2::RED) {
    //     const Subcurve* sc = *it;
    //     std::cout << "  Red: "
    //               << sc->red_halfedge_handle()->source()->point() << " => "
    //               << sc->red_halfedge_handle()->target()->point() << ", "
    //               << sc->red_halfedge_handle()->direction()
    //               << std::endl;
    //   }
    //   else {
    //     const Subcurve* sc = *it;
    //     std::cout << "  Blue: "
    //               << sc->blue_halfedge_handle()->source()->point() << " => "
    //               << sc->blue_halfedge_handle()->target()->point() << ", "
    //               << sc->blue_halfedge_handle()->direction()
    //               << std::endl;
    //   }
    // }

    if (event->parameter_space_in_y() == ARR_TOP_BOUNDARY) {
      // The curve is incident to the north pole; therefore,
      //  (i) the event has only left curves, and
      // (ii) the event point is the rightend (ARR_MIN_END) of all left curves.
      CGAL_assertion(event->number_of_right_curves() == 0);
      CGAL_assertion(event->number_of_left_curves() != 0);
      Subcurve_reverse_iterator it_end = event->left_curves_rend();

      // Handle red curves
      Subcurve_reverse_iterator it_red = event->left_curves_rbegin();
      while ((it_red != it_end) && ((*it_red)->color() == Gt2::BLUE)) ++it_red;
      if (it_red != it_end) {
        const Subcurve* sc = *it_red;
        CGAL_assertion(sc->color() != Gt2::BLUE);
        m_red_nf = sc->red_halfedge_handle()->face();
      }

      // Handle blue curves
      Subcurve_reverse_iterator it_blue = event->left_curves_rbegin();
      while ((it_blue != it_end) && ((*it_blue)->color() == Gt2::RED)) ++it_blue;
      if (it_blue != it_end) {
        const Subcurve* sc = *it_blue;
        CGAL_assertion(sc->color() != Gt2::RED);
        m_blue_nf = sc->blue_halfedge_handle()->face();
      }

      return;
    }

    // The curve is incident to, or lies on, the identification curve.
    Subcurve_reverse_iterator itr_end(event->right_curves_rend());

    // Handle red curves
    // First use the right curves if exists.
    Subcurve_reverse_iterator itr(event->right_curves_rbegin());
    while ((itr != itr_end) && ((*itr)->color() == Gt2::BLUE)) ++itr;
    if (itr != itr_end) {
      const Subcurve* sc = *itr;
      // Case 1.1.
      // The right event of a right curve coincides with the identification
      // curve. It implies that the entire curve lies on the identification
      // curve. In this case, the desired face is incident to the halfedge.
      //
      //             o : the right event
      //             |
      //             |_
      //             |/
      //  the event: o
      //             |
      //
      // Case 1.2.
      // The curve extends to the internal parameter space in X.
      // In this case, the desired face is incident to the twin halfedge.
      //
      //             |
      //             |
      //             |----->o : the right event
      //  the event: o<-----
      //             |
      //
      // Case 2.
      // The event does not have appropriately colored right-curves.
      // It has a left curve; it must lie on the identification curve.
      //
      //             o : the event
      //             |
      //             |_
      //             |/
      //             o
      //             |
      //
      auto* right_event = sc->right_event();
      m_red_nf = (right_event->parameter_space_in_x() == ARR_LEFT_BOUNDARY) ?
        sc->red_halfedge_handle()->face() :
        sc->red_halfedge_handle()->twin()->face();
    }
    else {
      Subcurve_reverse_iterator itl(event->left_curves_rbegin());
      Subcurve_reverse_iterator itl_end(event->left_curves_rend());
      while ((itl != itl_end) && ((*itl)->color() == Gt2::BLUE)) ++itl;
      if (itl != itl_end) {
        const Subcurve* sc = *itl;
        m_red_nf = sc->red_halfedge_handle()->face();
      }
    }

    // Handle blue curves
    // First use the right curves if exists.
    itr = event->right_curves_rbegin();
    while ((itr != itr_end) && ((*itr)->color() == Gt2::RED)) ++itr;
    if (itr != itr_end) {
      const Subcurve* sc = *itr;
      auto* right_event = sc->right_event();
      m_blue_nf = (right_event->parameter_space_in_x() == ARR_LEFT_BOUNDARY) ?
        sc->blue_halfedge_handle()->face() :
        sc->blue_halfedge_handle()->twin()->face();
    }
    else {
      Subcurve_reverse_iterator itl(event->left_curves_rbegin());
      Subcurve_reverse_iterator itl_end(event->left_curves_rend());
      while ((itl != itl_end) && ((*itl)->color() == Gt2::RED)) ++itl;
      if (itl != itl_end) {
        const Subcurve* sc = *itl;
        m_blue_nf = sc->blue_halfedge_handle()->face();
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
