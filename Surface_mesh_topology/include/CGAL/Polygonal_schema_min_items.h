// Copyright (c) 2019 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CGAL_POLYGONAL_SCHEMA_MIN_ITEMS_H
#define CGAL_POLYGONAL_SCHEMA_MIN_ITEMS_H 1

#include <CGAL/license/Surface_mesh_topology.h>
#include <string>

namespace CGAL {
namespace Surface_mesh_topology {

  /** @file Polygonal_schema_min_items.h
   * Definition of min item class for Polygonal_schema.
   */

  /** Minimal items for polygonal schema.
   * Generic_map_min_items defines what is the minimal item class for a generic map.
   * One struct associated with darts, having one std::string named m_label..
   */
  struct Polygonal_schema_min_items
  {
    template < class Refs >
    struct Dart_wrapper
    {
      struct Info_for_darts
      { std::string m_label; };

      typedef Info_for_darts Dart_info;
    };
  };

} // namespace Surface_mesh_topology
} // namespace CGAL

#endif // CGAL_POLYGONAL_SCHEMA_MIN_ITEMS_H
// EOF //
