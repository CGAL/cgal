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
// Author(s)     : Tali Zvi <talizvi@post.tau.ac.il>
//                 Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_ARR_OVERLAY_EVENT_H
#define CGAL_ARR_OVERLAY_EVENT_H

#include <CGAL/license/Arrangement_on_surface_2.h>

/*! \file
 *
 * Definition of the Arr_construction_event_base class-template.
 */

#include <CGAL/Surface_sweep_2/Arr_construction_event_base.h>
#include <CGAL/Surface_sweep_2/Arr_overlay_subcurve.h>

namespace CGAL {

/* \class Arr_overlay_event
 *
 * This template represents an event used by the surface-sweep framework, where
 * the input curves for the surface-sweep procedure are extracted from two
 * arrangements that are overlaid.
 *
 * \tparam GeometryTraits_2 the geometry traits.
 * \tparam Arrangement_ the type of the arrangement that is the resulting
 *                      arrangement the overlay process.
 * \tparam Allocator_ a type of an element that is used to acquire/release
 *                    memory for elements of the event queue and the status
 *                    structure, and to construct/destroy the elements in that
 *                    memory. The type must meet the requirements of Allocator.
 */
template <typename GeometryTraits_2, typename Arrangement_,
          typename Allocator_ = CGAL_ALLOCATOR(int)>
class Arr_overlay_event :
  public Arr_construction_event_base<GeometryTraits_2,
                                     Arr_overlay_subcurve<GeometryTraits_2,
                                                          Arr_overlay_event<
                                                            GeometryTraits_2,
                                                            Arrangement_,
                                                            Allocator_>,
                                                          Allocator_>,
                                     Arrangement_>
{
public:
  /*! Construct default. */
  Arr_overlay_event() {}
};

} // namespace CGAL

#endif
