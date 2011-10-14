// Copyright (c) 2010 CNRS, LIRIS, http://liris.cnrs.fr/, All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
#ifndef CGAL_COMBINATORIAL_MAP_WITH_POINTS_MIN_ITEMS_H
#define CGAL_COMBINATORIAL_MAP_WITH_POINTS_MIN_ITEMS_H 1

#include <CGAL/Dart.h>
#include <CGAL/Cell_attribute_with_point.h>

namespace CGAL {

/** @file Combinatorial_map_with_points_min_items.h
 * Definition of min item class for map with points.
 */

/** Minimal items for combinatorial map with points.
 * Combinatorial_map_with_points_min_items defines what is the item class
 * for a map with points. It provides definitions for attributes associated
 * to vertices (containing points), and darts. The traits class must provide the
 * respective type for the point.
 */
  template <unsigned int d1
						/*,
												unsigned int d2=d1, 
													class Traits_=typename Default_template_argument<d2>::type*/ >
  struct Combinatorial_map_with_points_min_items
  {
		//    typedef Traits_           Traits;             

    /// Dart_wrapper defines the type of darts used.
    template <class Refs>
    struct Dart_wrapper
    {
      typedef CGAL::Dart<d1, Refs> Dart;

      typedef Cell_attribute_with_point<Refs> Vertex_attrib;    
      typedef CGAL::cpp0x::tuple<Vertex_attrib> Attributes;
    };
  };

} // namespace CGAL

#endif // CGAL_COMBINATORIAL_MAP_WITH_POINTS_MIN_ITEMS_H //
// EOF //
