// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Tali Zvi <talizvi@post.tau.ac.il>,
//             Baruch Zukerman <baruchzu@post.tau.ac.il>
//             Ron Wein <wein@post.tau.ac.il>
//             Efi Fogel <efif@gmail.com>

#ifndef CGAL_SURFACE_SWEEP_2_DEFAULT_EVENT_H
#define CGAL_SURFACE_SWEEP_2_DEFAULT_EVENT_H

#include <CGAL/license/Surface_sweep_2.h>

/*! \file
 *
 * Defintion of the Default_event class.
 */

#include <CGAL/Surface_sweep_2/Default_event_base.h>
#include <CGAL/Surface_sweep_2/Default_subcurve.h>

namespace CGAL {
namespace Surface_sweep_2 {

/*! \class Default_event
 *
 * This template represents an event used by the surface-sweep framework, where
 * the input curves for the surface-sweep procedure may overlap.
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
 * Default_event_base; do not use this class as base in your derivation.
 */
template <typename GeometryTraits_2,
          typename Allocator_ = CGAL_ALLOCATOR(int)>
class Default_event :
  public Default_event_base<GeometryTraits_2,
                            Default_subcurve<GeometryTraits_2,
                                             Default_event<GeometryTraits_2,
                                                           Allocator_>,
                                             Allocator_> >
{
public:
  /*! Construct default. */
  Default_event() {}
};

} // namespace Surface_sweep_2
} // namespace CGAL

#endif
