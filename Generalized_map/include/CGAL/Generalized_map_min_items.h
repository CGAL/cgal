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
#ifndef CGAL_GENERALIZED_MAP_MIN_ITEMS_H
#define CGAL_GENERALIZED_MAP_MIN_ITEMS_H 1

#include <CGAL/Generalized_map_dart.h>

namespace CGAL {

/** @file Generalized_map_min_items.h
 * Definition of min item class for dD generalized map.
 */

/** Minimal items for dD generalized map.
 * Generalized_map_min_items defines what is the minimal item
 * class for a d-gmap It provides definitions for darts without attribute.
 */
template <unsigned int d>
struct Generalized_map_min_items
{
   /// Dart_wrapper defines the type of darts used, and enabled attributes.
   template < class Refs >
   struct Dart_wrapper
   {
     typedef CGAL::Generalized_map_dart< d, Refs > Dart;
     typedef CGAL::cpp11::tuple<> Attributes;
   };
};

} // namespace CGAL

#endif // CGAL_GENERALIZED_MAP_MIN_ITEMS_H //
// EOF //
