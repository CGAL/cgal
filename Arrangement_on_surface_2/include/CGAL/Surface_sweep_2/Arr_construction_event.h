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

template <typename GeometryTraits_2, typename Arrangement_,
          template <typename, typename>
          class SurfaceSweepBaseEvent = Ss2::Default_event_base,
          template <typename, typename, typename>
          class SurfaceSweepBaseCurve = Ss2::Default_subcurve>
class Arr_construction_event :
  public Arr_construction_event_base<
    GeometryTraits_2,
    Arr_construction_subcurve<GeometryTraits_2,
                              Arr_construction_event<GeometryTraits_2,
                                                     Arrangement_,
                                                     SurfaceSweepBaseEvent,
                                                     SurfaceSweepBaseCurve>,
                              SurfaceSweepBaseCurve>,
    Arrangement_, SurfaceSweepBaseEvent>
{
public:
  /*! Construct default. */
  Arr_construction_event() {}
};

} // namespace CGAL

#endif
