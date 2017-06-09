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
 */
template <typename GeometryTraits_2>
class Default_event :
  public Default_event_base<GeometryTraits_2,
                            Default_subcurve<GeometryTraits_2,
                                             Default_event<GeometryTraits_2> > >
{
public:
  /*! Construct default. */
  Default_event() {}
};

} // namespace Surface_sweep_2
} // namespace CGAL

#endif
