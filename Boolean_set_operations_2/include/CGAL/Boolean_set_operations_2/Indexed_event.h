// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// Author(s) : Baruch Zukerman <baruchzu@post.tau.ac.il>
//             Ron Wein <wein@post.tau.ac.il>
//             Efi Fogel <efifogel@gmail.com>

#ifndef CGAL_BSO_2_INDEXED_VISITOR_H
#define CGAL_BSO_2_INDEXED_VISITOR_H

#include <CGAL/license/Boolean_set_operations_2.h>

#include <CGAL/Surface_sweep_2/Arr_construction_event_base.h>

namespace CGAL {

/* \class Indexed_event
 */
template <typename GeometryTraits_2, typename Arrangement_,
          typename Allocator_ = CGAL_ALLOCATOR(int)>
class Indexed_event :
  public Arr_construction_event_base<
    GeometryTraits_2,
    Arr_construction_subcurve<GeometryTraits_2,
                              Indexed_event<GeometryTraits_2,
                                            Arrangement_,
                                            Allocator_>,
                              Allocator_>,
    Arrangement_>
{
private:
  unsigned int m_index;

public:
  Indexed_event() : m_index (0) {}

  unsigned int index() const { return (m_index); }

  void set_index(unsigned int index) { m_index = index; }
};

} // namespace CGAL

#endif
