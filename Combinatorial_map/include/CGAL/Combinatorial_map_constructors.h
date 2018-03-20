// Copyright (c) 2010-2011 CNRS and LIRIS' Establishments (France).
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
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CGAL_COMBINATORIAL_MAP_CONSTRUCTORS_H
#define CGAL_COMBINATORIAL_MAP_CONSTRUCTORS_H 1

#include <CGAL/config.h>

#ifndef CGAL_NO_DEPRECATED_CODE

namespace CGAL {

  /** @file Combinatorial_map_constructors.h
   * Basic creation operations  for a combinatorial map.
   * Create edge, triangle, quadrilateral, tetrahedron, hexahedron.
   */

  /** Create an edge.
   * @param amap the used combinatorial map.
   * @return a dart of the new edge.
   */
  template < class Map >
  CGAL_DEPRECATED typename Map::Dart_handle make_edge(Map& amap)
  { return amap.make_edge(); }

  /** Create a combinatorial polygon of length alg
   * (a cycle of alg darts beta1 links together).
   * @param amap the used combinatorial map.
   * @return a new dart.
   */
  template < class Map >
  CGAL_DEPRECATED typename Map::Dart_handle make_combinatorial_polygon(Map& amap,
                                                       unsigned int alg)
  { return amap.make_combinatorial_polygon(alg); }

  /** Test if a face is a combinatorial polygon of length alg
   *  (a cycle of alg darts beta1 links together).
   * @param amap the used combinatorial map.
   * @param adart an intial dart
   * @return true iff the face containing adart is a polygon of length alg.
   */
  template < class Map >
  CGAL_DEPRECATED bool is_face_combinatorial_polygon(const Map& amap,
                                     typename Map::Dart_const_handle adart,
                                     unsigned int alg)
  { return amap.is_face_combinatorial_polygon(adart, alg); }

  /** Create a combinatorial tetrahedron from 4 triangles.
   * @param amap the used combinatorial map.
   * @param d1 a dart onto a first triangle.
   * @param d2 a dart onto a second triangle.
   * @param d3 a dart onto a third triangle.
   * @param d4 a dart onto a fourth triangle.
   * @return a new dart.
   */
  template < class Map >
  typename Map::Dart_handle
  CGAL_DEPRECATED make_combinatorial_tetrahedron(Map& amap,
                                 typename Map::Dart_handle d1,
                                 typename Map::Dart_handle d2,
                                 typename Map::Dart_handle d3,
                                 typename Map::Dart_handle d4)
  { return amap.make_combinatorial_tetrahedron(d1, d2, d3, d4); }

  /** Test if a volume is a combinatorial tetrahedron.
   * @param amap the used combinatorial map.
   * @param adart an intial dart
   * @return true iff the volume containing adart is a combinatorial tetrahedron.
   */
  template < class Map >
  CGAL_DEPRECATED bool is_volume_combinatorial_tetrahedron(const Map& amap,
                                           typename Map::Dart_const_handle d1)
  { return amap.is_volume_combinatorial_tetrahedron(d1); }

  /** Create a new combinatorial tetrahedron.
   * @param amap the used combinatorial map.
   * @return a new dart.
   */
  template < class Map >
  CGAL_DEPRECATED typename Map::Dart_handle make_combinatorial_tetrahedron(Map& amap)
  { return amap.make_combinatorial_tetrahedron(); }

  /** Create a combinatorial hexahedron from 6 quadrilaterals.
   * @param amap the used combinatorial map.
   * @param d1 a dart onto a first quadrilateral.
   * @param d2 a dart onto a second quadrilateral.
   * @param d3 a dart onto a third quadrilateral.
   * @param d4 a dart onto a fourth quadrilateral.
   * @param d5 a dart onto a fifth quadrilateral.
   * @param d6 a dart onto a sixth quadrilateral.
   * @return a dart of the new cuboidal_cell.
   */
  template < class Map >
  typename Map::Dart_handle
  CGAL_DEPRECATED make_combinatorial_hexahedron(Map& amap,
                                typename Map::Dart_handle d1,
                                typename Map::Dart_handle d2,
                                typename Map::Dart_handle d3,
                                typename Map::Dart_handle d4,
                                typename Map::Dart_handle d5,
                                typename Map::Dart_handle d6)
  { return amap.make_combinatorial_hexahedron(d1, d2, d3, d4, d5, d6); }

  /** Test if a volume is a combinatorial hexahedron.
   * @param amap the used combinatorial map.
   * @param adart an intial dart
   * @return true iff the volume containing adart is a combinatorial hexahedron.
   */
  template < class Map >
  CGAL_DEPRECATED bool is_volume_combinatorial_hexahedron(const Map& amap,
                                          typename Map::Dart_const_handle d1)
  { return amap.is_volume_combinatorial_hexahedron(d1); }

  /** Create a new combinatorial hexahedron.
   * @param amap the used combinatorial map.
   * @return a new dart.
   */
  template < class Map >
  CGAL_DEPRECATED typename Map::Dart_handle make_combinatorial_hexahedron(Map& amap)
  { return amap.make_combinatorial_hexahedron(); }

} // namespace CGAL

#endif // CGAL_NO_DEPRECATED_CODE

#endif // CGAL_COMBINATORIAL_MAP_CONSTRUCTORS_H //
// EOF //
