// Copyright (c) 2007,2009,2010,2011,2013 Tel-Aviv University (Israel).
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
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Ron Wein <wein@post.tau.ac.il>
//                 Efi Fogel <efif@post.tau.ac.il>
//

#ifndef CGAL_ARR_SPHERICAL_BATCHED_PL_HELPER_H
#define CGAL_ARR_SPHERICAL_BATCHED_PL_HELPER_H

#include <CGAL/license/Arrangement_on_surface_2.h>


/*! \file
 * Definition of the Arr_spherical_batched_pl_helper class-template.
 */

namespace CGAL {

#include <CGAL/Sweep_line_empty_visitor.h>

/*! \class Arr_spherical_batched_pl_helper
 * A helper class for the batched point-location sweep-line visitor, suitable
 * for an Arrangement_on_surface_2 instantiated with a topology-traits class
 * for bounded curves in the plane.
 */
template <class Traits_, class Arrangement_>
class Arr_spherical_batched_pl_helper
{
public:

  typedef Traits_                                      Traits_2;
  typedef Arrangement_                                 Arrangement_2;

  typedef typename Arrangement_2::Face_const_handle    Face_const_handle;

  typedef Sweep_line_empty_visitor<Traits_2>           Base_visitor;
  typedef typename Base_visitor::Event                 Event;
  typedef typename Base_visitor::Subcurve              Subcurve;
  typedef typename Event::Subcurve_iterator            Subcurve_iterator;

protected:

  typedef typename Arrangement_2::Topology_traits      Topology_traits;

  // Data members:
  //! The topology-traits class.
  const Topology_traits* m_top_traits;

  //! The unbounded arrangement face.
  Face_const_handle m_spherical_face;

public:
  /*! Constructor.
   * \param arr The arrangement.
   */
  Arr_spherical_batched_pl_helper(const Arrangement_2* arr) :
    m_top_traits(arr->topology_traits())
  {}

  /// \name Notification functions.
  //@{

  /*! A notification issued before the sweep process starts. */
  void before_sweep()
  { m_spherical_face = Face_const_handle(m_top_traits->spherical_face()); }

  /*! A notification invoked after the sweep-line finishes handling the given
   * event.
   */
  void after_handle_event(Event* event) {

    if (event->parameter_space_in_y() == ARR_TOP_BOUNDARY) {
      Arr_curve_end ind = ((event->number_of_left_curves() == 0) &&
                           (event->number_of_right_curves() != 0)) ?
        ARR_MIN_END : ARR_MAX_END;
      Subcurve_iterator it, nit, it_end;
      if (ind == ARR_MIN_END) {
        it = nit = event->right_curves_begin();
        it_end = event->right_curves_end();
      } else {
        it = nit = event->left_curves_begin();
        it_end = event->left_curves_end();
      }

      ++nit;
      if (it != it_end) {
        while (nit != it_end) {
          ++it;
          ++nit;
        }
      }
      const Subcurve* sc = *it;
      // pick the one facing the top right corner now
      m_spherical_face = sc->last_curve().halfedge_handle()->face();
    }

    // EBEB 2013-12-012 do similar stuff for right boundary
    if (event->parameter_space_in_y() == ARR_RIGHT_BOUNDARY) {
      Subcurve_iterator it, nit, it_end;
      it = nit = event->left_curves_begin();
      it_end = event->left_curves_end();

      ++nit;
      if (it != it_end) {
        while (nit != it_end) {
          ++it;
          ++nit;
        }
      }
      const Subcurve* sc = *it;
      // pick the one facing the top right corner now
      CGAL_assertion(sc->last_curve().halfedge_handle()->direction() == ARR_LEFT_TO_RIGHT);
      m_spherical_face = sc->last_curve().halfedge_handle()->face();
    }


    return;
  }
  //@}

  /*! Obtain the current top face. */
  Face_const_handle top_face() const { return m_spherical_face; }
};

} //namespace CGAL

#endif
