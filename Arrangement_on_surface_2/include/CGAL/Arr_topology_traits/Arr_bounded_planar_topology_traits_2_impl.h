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
// Author(s)     : Ron Wein <wein@post.tau.ac.il>
//                 Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_ARR_BOUNDED_PLANAR_TOPOLOGY_TRAITS_2_IMPL_H
#define CGAL_ARR_BOUNDED_PLANAR_TOPOLOGY_TRAITS_2_IMPL_H

#include <CGAL/license/Arrangement_on_surface_2.h>

/*! \file
 *
 * Member-function definitions for the
 * Arr_bounded_planar_topology_traits_2<GeomTraits> class.
 */

namespace CGAL {

//-----------------------------------------------------------------------------
// Assign the contents of another topology-traits class.
//
template <typename GeometryTraits_2, typename Dcel_>
void Arr_bounded_planar_topology_traits_2<GeometryTraits_2, Dcel_>::
assign(const Self& other)
{
  // Assign the base class.
  Base::assign(other);

  // Update the topology-traits properties after the DCEL have been updated.
  dcel_updated();
}

//-----------------------------------------------------------------------------
// Initialize an empty DCEL structure.
//
template <typename GeometryTraits_2, typename Dcel_>
void Arr_bounded_planar_topology_traits_2<GeometryTraits_2, Dcel_>::init_dcel()
{
  // Clear the current DCEL.
  this->m_dcel.delete_all();

  // Create the unbounded face.
  unb_face = this->m_dcel.new_face();

  unb_face->set_unbounded(true);
  unb_face->set_fictitious(false);
}

//-----------------------------------------------------------------------------
// Make the necessary updates after the DCEL structure have been updated.
//
template <typename GeometryTraits_2, typename Dcel_>
void Arr_bounded_planar_topology_traits_2<GeometryTraits_2, Dcel_>::
dcel_updated()
{
  // Go over the DCEL faces and locate the unbounded face.
  unb_face = nullptr;
  typename Dcel::Face_iterator fit = this->m_dcel.faces_begin();
  for (; fit != this->m_dcel.faces_end(); ++fit) {
    if (fit->is_unbounded()) {
      unb_face = &(*fit);
      break;
    }
  }
  CGAL_assertion(unb_face != nullptr);
}

} // namespace CGAL

#endif
