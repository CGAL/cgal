// Copyright (c) 2019 CNRS and LIRIS' Establishments (France).
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
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CGAL_POLYGONAL_SCHEMA_MIN_ITEMS_H
#define CGAL_POLYGONAL_SCHEMA_MIN_ITEMS_H 1

#include <CGAL/license/Surface_mesh_topology.h>

#include <cstring>
  
namespace CGAL {
namespace Surface_mesh_topology {

  /** @file Polygonal_schema_min_items.h
   * Definition of min item class for Polygonal_schema.
   */

  /** Minimal items for polygonal schema.
   * Generic_map_min_items defines what is the minimal item class for a generic map.
   * One struct associated with darts, having one char* named m_label..
   */
  struct Polygonal_schema_min_items
  {
    template < class Refs >
    struct Dart_wrapper
    {
      struct Info_for_darts
      {
        char* m_label;

        Info_for_darts() : m_label(nullptr)
        {}
      };
      
      typedef Info_for_darts Dart_info;
    };
  };

} // namespace Surface_mesh_topology
} // namespace CGAL

#endif // CGAL_POLYGONAL_SCHEMA_MIN_ITEMS_H
// EOF //
