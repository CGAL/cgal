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
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s) : Tali Zvi <talizvi@post.tau.ac.il>,
//             Baruch Zukerman <baruchzu@post.tau.ac.il>
//             Ron Wein <wein@post.tau.ac.il>
//             Efi Fogel <efifogel@gmail.com>

#ifndef CGAL_SURFACE_SWEEP_2_NO_OVERLAP_EVENT_H
#define CGAL_SURFACE_SWEEP_2_NO_OVERLAP_EVENT_H

#include <CGAL/license/Surface_sweep_2.h>

/*! \file
 *
 * Defintion of the No_overlap_event class.
 */

#include <CGAL/Surface_sweep_2/No_overlap_event_base.h>
#include <CGAL/Surface_sweep_2/No_overlap_subcurve.h>

namespace CGAL {
namespace Surface_sweep_2 {

/*! \class No_overlap_event
 *
 * This template represents an event used by the surface-sweep framework, where
 * the input curves for the surface-sweep procedure are guaranteed not to
 * overlap.
 *
 * \tparam GeometryTraits_2 the geometry traits.
 * \tparam Allocator_ a type of an element that is used to acquire/release
 *                    memory for elements of the event queue and the status
 *                    structure, and to construct/destroy the elements in that
 *                    memory. The type must meet the requirements of Allocator.
 *
 * We exploit the curiously recurring template pattern (CRTP) idiom to establish
 * an interdependency between the curve and the event types, which are template
 * parameters of the surface-sweep visitor class templates. It enables the
 * definition of these two types, which refer one to another; (the curves to the
 * right of an event and the curves to its left are data members of the event,
 * and the two events associated with the endpoints of a curve are data memebrs
 * of the curve.)
 *
 * If you need to represent an event with additional data members, introduce a
 * new type, say x, that derives from x_base, and have x_base derive from
 * No_overlap_event_base; do not use this class as base in your derivation.
 */
template <typename GeometryTraits_2,
          typename Allocator_ = CGAL_ALLOCATOR(int)>
class No_overlap_event :
  public No_overlap_event_base<GeometryTraits_2,
                               No_overlap_subcurve<GeometryTraits_2,
                                                   No_overlap_event<
                                                     GeometryTraits_2,
                                                     Allocator_>,
                                                   Allocator_> >
{
public:
  /*! Construct default. */
  No_overlap_event() {}
};

} // namespace Surface_sweep_2
} // namespace CGAL

#endif
