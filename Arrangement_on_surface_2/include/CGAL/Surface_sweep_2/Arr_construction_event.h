// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Tali Zvi <talizvi@post.tau.ac.il>
//                 Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_ARR_CONSTRUCTION_EVENT_H
#define CGAL_ARR_CONSTRUCTION_EVENT_H

#include <CGAL/license/Arrangement_on_surface_2.h>

/*! \file
 *
 * Definition of the Arr_construction_event class-template.
 */

#include <CGAL/Surface_sweep_2/Default_event_base.h>
#include <CGAL/Surface_sweep_2/Arr_construction_event_base.h>
#include <CGAL/Surface_sweep_2/Arr_construction_subcurve.h>
#include <CGAL/Surface_sweep_2/Default_subcurve.h>

namespace CGAL {

namespace Ss2 = Surface_sweep_2;

/*! \class Default_event
 *
 * This template represents an event used by the surface-sweep framework, where
 * the input curves for the surface-sweep procedure are to be inserted into a
 * 2D arrangement.
 *
 * \tparam GeometryTraits_2 the geometry traits.
 * \tparam Arrangement_ the type of the costructed arrangement.
 * \tparam Allocator_ a type of an element that is used to acquire/release
 *                    memory for elements of the event queue and the status
 *                    structure, and to construct/destroy the elements in that
 *                    memory. The type must meet the requirements of Allocator.
 * \tparam SurfaceSweepBaseEvent the base class of the event.
 * \tparam SurfaceSweepBaseCurve the base class of the subcurve.
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
 * Arr_construction_event_base; do not use this class as base in your
 * derivation.
 */
template <typename GeometryTraits_2, typename Arrangement_,
          typename Allocator_ = CGAL_ALLOCATOR(int),
          template <typename, typename>
          class SurfaceSweepBaseEvent = Ss2::Default_event_base,
          template <typename, typename, typename, typename>
          class SurfaceSweepBaseCurve = Ss2::Default_subcurve>
class Arr_construction_event :
  public Arr_construction_event_base<
    GeometryTraits_2,
    Arr_construction_subcurve<GeometryTraits_2,
                              Arr_construction_event<GeometryTraits_2,
                                                     Arrangement_,
                                                     Allocator_,
                                                     SurfaceSweepBaseEvent,
                                                     SurfaceSweepBaseCurve>,
                              Allocator_,
                              SurfaceSweepBaseCurve>,
    Arrangement_, SurfaceSweepBaseEvent>
{
public:
  /*! Construct default. */
  Arr_construction_event() {}
};

} // namespace CGAL

#endif
