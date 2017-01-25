// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
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
// Author(s)     : Ron Wein <wein@post.tau.ac.il>
//                 Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_ARR_BOUNDED_PLANAR_TOPOLOGY_TRAITS_2_IMPL_H
#define CGAL_ARR_BOUNDED_PLANAR_TOPOLOGY_TRAITS_2_IMPL_H

#include <CGAL/license/Arrangement_on_surface_2.h>


/*! \file
 * Member-function definitions for the
 * Arr_bounded_planar_topology_traits_2<GeomTraits> class.
 */

namespace CGAL {

//-----------------------------------------------------------------------------
// Assign the contents of another topology-traits class.
//
template <typename GeomTraits_, typename Dcel_>
void Arr_bounded_planar_topology_traits_2<GeomTraits_, Dcel_>::
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
template <typename GeomTraits_, typename Dcel_>
void Arr_bounded_planar_topology_traits_2<GeomTraits_, Dcel_>::init_dcel()
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
template <typename GeomTraits_, typename Dcel_>
void Arr_bounded_planar_topology_traits_2<GeomTraits_, Dcel_>::dcel_updated()
{
  // Go over the DCEL faces and locate the unbounded face.
  unb_face = NULL;
  typename Dcel::Face_iterator fit = this->m_dcel.faces_begin();
  for (; fit != this->m_dcel.faces_end(); ++fit) {
    if (fit->is_unbounded()) {
      unb_face = &(*fit);
      break;
    }
  }
  CGAL_assertion(unb_face != NULL);
}

} //namespace CGAL

#endif
