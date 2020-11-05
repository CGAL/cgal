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

#ifndef CGAL_ARR_BOUNDED_PLANAR_BATCHED_PL_HELPER_H
#define CGAL_ARR_BOUNDED_PLANAR_BATCHED_PL_HELPER_H

#include <CGAL/license/Arrangement_on_surface_2.h>

/*! \file
 *
 * Definition of the Arr_bounded_planar_batched_pl_helper class-template.
 */

namespace CGAL {

/*! \class Arr_bounded_planar_batched_pl_helper
 *
 * A helper class for the batched point-location surface-sweep visitor, suitable
 * for an Arrangement_on_surface_2 instantiated with a topology-traits class
 * for bounded curves in the plane.
 */
template <typename GeometryTraits_2, typename Arrangement_, typename Event_,
          typename Subcurve_>
class Arr_bounded_planar_batched_pl_helper {
public:
  typedef GeometryTraits_2                              Geometry_traits_2;
  typedef Arrangement_                                  Arrangement_2;
  typedef Event_                                        Event;
  typedef Subcurve_                                     Subcurve;
  typedef typename Subcurve::Allocator                  Allocator;

private:
  typedef Geometry_traits_2                             Gt2;
  typedef typename Arrangement_2::Face_const_handle     Face_const_handle;

public:
  typedef typename Arrangement_2::Topology_traits       Topology_traits;

protected:
  // Data members:
  const Topology_traits* m_top_traits;  // The topology-traits class.
  Face_const_handle m_unb_face;         // The unbounded arrangement face.

public:
  /*! Constructor.
   * \param arr The arrangement.
   */
  Arr_bounded_planar_batched_pl_helper(const Arrangement_2* arr) :
    m_top_traits(arr->topology_traits())
  {}

  /// \name Notification functions.
  //@{

  /* A notification issued before the sweep process starts. */
  void before_sweep()
  {
    // Get the unbounded face.
    m_unb_face = Face_const_handle(m_top_traits->unbounded_face());
  }

  /*! A notification invoked after the sweep-line finishes handling the given
   * event.
   */
  void after_handle_event(Event* /* event */) { return; }
  //@}

  /*! Get the current top face. */
  Face_const_handle top_face() const { return (m_unb_face); }
};

} // namespace CGAL

#endif
