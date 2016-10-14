// Copyright (c) 2016 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CGAL_GENERIC_MAP_MIN_ITEMS_H
#define CGAL_GENERIC_MAP_MIN_ITEMS_H 1

namespace CGAL {

  /** @file Generic_map_min_items.h
   * Definition of min item class for dD generic map.
   */

  /** Minimal items for dD generic map.
   * Generic_map_min_items defines what is the minimal item class for a generic map.
   * No information associated with darts, no enabled attribute.
   */
  struct Generic_map_min_items
  {
    template < class Refs >
    struct Dart_wrapper
    {
      typedef void Dart_info;
      typedef CGAL::cpp11::tuple<> Attributes;
    };
  };

} // namespace CGAL

#endif // CGAL_GENERIC_MAP_MIN_ITEMS_H
// EOF //
