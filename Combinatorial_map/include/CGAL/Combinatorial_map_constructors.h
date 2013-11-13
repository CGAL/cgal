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
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CGAL_COMBINATORIAL_MAP_CONSTRUCTORS_H
#define CGAL_COMBINATORIAL_MAP_CONSTRUCTORS_H 1

#include <CGAL/Combinatorial_map_basic_operations.h>

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
  typename Map::Dart_handle make_edge(Map& amap)
  {
    typename Map::Dart_handle d1 = amap.create_dart();
    typename Map::Dart_handle d2 = amap.create_dart();
    amap.basic_link_beta_for_involution(d1, d2, 2);
    return d1;
  }

  /** Create a combinatorial polygon of length alg 
   * (a cycle of alg darts beta1 links together).
   * @param amap the used combinatorial map.
   * @return a new dart.
   */
  template < class Map >
  typename Map::Dart_handle make_combinatorial_polygon(Map& amap, 
                                                       unsigned int alg)
  {
    CGAL_assertion(alg>0);
  
    typename Map::Dart_handle start = amap.create_dart();
    typename Map::Dart_handle prev = start;
    for ( unsigned int nb=1; nb<alg; ++nb )
    {
      typename Map::Dart_handle cur = amap.create_dart();
      amap.basic_link_beta_1(prev, cur);
      prev=cur;
    }
  
    amap.basic_link_beta_1(prev, start);
    return start;
  }

/** Test if a face is a combinatorial polygon of length alg
 *  (a cycle of alg darts beta1 links together).
 * @param amap the used combinatorial map.
 * @param adart an intial dart
 * @return true iff the face containing adart is a polygon of length alg.
 */
template < class Map >
bool is_face_combinatorial_polygon(const Map& amap,
                                   typename Map::Dart_const_handle adart,
                                   unsigned int alg)
{
  CGAL_assertion(alg>0);

  unsigned int nb = 0;
  typename Map::Dart_const_handle cur = adart;
  do
  {
    ++nb;
    if ( cur==amap.null_dart_handle ) return false; // Open face
    cur = amap.beta(cur,1);
  }
  while( cur!=adart );
  return (nb==alg);
}

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
  make_combinatorial_tetrahedron(Map& amap,
                                 typename Map::Dart_handle d1,
                                 typename Map::Dart_handle d2,
                                 typename Map::Dart_handle d3,
                                 typename Map::Dart_handle d4)
  {
    amap.basic_link_beta_for_involution(d1, d2, 2);
    amap.basic_link_beta_for_involution(d3, amap.beta(d2, 0), 2);
    amap.basic_link_beta_for_involution(amap.beta(d1, 1), amap.beta(d3, 0), 2);
    amap.basic_link_beta_for_involution(d4, amap.beta(d2, 1), 2);
    amap.basic_link_beta_for_involution(amap.beta(d4, 0), amap.beta(d3, 1), 2);
    amap.basic_link_beta_for_involution(amap.beta(d4, 1), amap.beta(d1, 0), 2);

    return d1;
  }

/** Test if a volume is a combinatorial tetrahedron.
 * @param amap the used combinatorial map.
 * @param adart an intial dart
 * @return true iff the volume containing adart is a combinatorial tetrahedron.
 */
template < class Map >
bool is_volume_combinatorial_tetrahedron(const Map& amap,
                                         typename Map::Dart_const_handle d1)
{
  typename Map::Dart_const_handle d2 = amap.beta(d1, 2);
  typename Map::Dart_const_handle d3 = amap.beta(d2, 0, 2);
  typename Map::Dart_const_handle d4 = amap.beta(d2, 1, 2);

  if ( d1==amap.null_dart_handle || d2==amap.null_dart_handle ||
       d3==amap.null_dart_handle || d4==amap.null_dart_handle ) return false;

  if ( !is_face_combinatorial_polygon(amap, d1, 3) ||
       !is_face_combinatorial_polygon(amap, d2, 3) ||
       !is_face_combinatorial_polygon(amap, d3, 3) ||
       !is_face_combinatorial_polygon(amap, d4, 3) ) return false;

  // TODO do better with marks (?).
  if ( belong_to_same_cell<Map,2,1>(amap, d1, d2) ||
       belong_to_same_cell<Map,2,1>(amap, d1, d3) ||
       belong_to_same_cell<Map,2,1>(amap, d1, d4) ||
       belong_to_same_cell<Map,2,1>(amap, d2, d3) ||
       belong_to_same_cell<Map,2,1>(amap, d2, d4) ||
       belong_to_same_cell<Map,2,1>(amap, d3, d4) ) return false;

  if ( amap.beta(d1,1,2)!=amap.beta(d3,0) ||
       amap.beta(d4,0,2)!=amap.beta(d3,1) ||
       amap.beta(d4,1,2)!=amap.beta(d1,0) ) return false;

  return true;
}


  /** Create a new combinatorial tetrahedron.
   * @param amap the used combinatorial map.
   * @return a new dart.
   */
  template < class Map >
  typename Map::Dart_handle make_combinatorial_tetrahedron(Map& amap)
  {
    typename Map::Dart_handle d1 = make_combinatorial_polygon(amap,3);
    typename Map::Dart_handle d2 = make_combinatorial_polygon(amap,3);
    typename Map::Dart_handle d3 = make_combinatorial_polygon(amap,3);
    typename Map::Dart_handle d4 = make_combinatorial_polygon(amap,3);

    return make_combinatorial_tetrahedron(amap, d1, d2, d3, d4);
  }

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
  make_combinatorial_hexahedron(Map& amap,
                                typename Map::Dart_handle d1,
                                typename Map::Dart_handle d2,
                                typename Map::Dart_handle d3,
                                typename Map::Dart_handle d4,
                                typename Map::Dart_handle d5,
                                typename Map::Dart_handle d6)
  {
    amap.basic_link_beta_for_involution(d1,
                                        amap.beta(d4, 1, 1), 2);
    amap.basic_link_beta_for_involution(amap.beta(d1, 1),
                                        amap.beta(d6, 0)    , 2);
    amap.basic_link_beta_for_involution(amap.beta(d1, 1, 1),
                                        d2                  , 2);
    amap.basic_link_beta_for_involution(amap.beta(d1, 0),
                                        d5                  , 2);

    amap.basic_link_beta_for_involution(d3,
                                        amap.beta(d2, 1, 1), 2);
    amap.basic_link_beta_for_involution(amap.beta(d3, 1),
                                        amap.beta(d6, 1)         , 2);
    amap.basic_link_beta_for_involution(amap.beta(d3, 1, 1),
                                        d4                  , 2);
    amap.basic_link_beta_for_involution(amap.beta(d3, 0),
                                        amap.beta(d5, 1, 1), 2);

    amap.basic_link_beta_for_involution(d6,
                                        amap.beta(d4, 1)         , 2);
    amap.basic_link_beta_for_involution(amap.beta(d6, 1, 1),
                                        amap.beta(d2, 1)         , 2);

    amap.basic_link_beta_for_involution(amap.beta(d5, 0),
                                        amap.beta(d4, 0)         , 2);
    amap.basic_link_beta_for_involution(amap.beta(d5, 1),
                                        amap.beta(d2, 0)         , 2);

    return d1;
  }

/** Test if a volume is a combinatorial hexahedron.
 * @param amap the used combinatorial map.
 * @param adart an intial dart
 * @return true iff the volume containing adart is a combinatorial hexahedron.
 */
template < class Map >
bool is_volume_combinatorial_hexahedron(const Map& amap,
                                        typename Map::Dart_const_handle d1)
{
  typename Map::Dart_const_handle d2 = amap.beta(d1, 1, 1, 2);
  typename Map::Dart_const_handle d3 = amap.beta(d2, 1, 1, 2);
  typename Map::Dart_const_handle d4 = amap.beta(d3, 1, 1, 2);
  typename Map::Dart_const_handle d5 = amap.beta(d1, 0, 2);
  typename Map::Dart_const_handle d6 = amap.beta(d4, 1, 2);

  if ( d1==amap.null_dart_handle || d2==amap.null_dart_handle ||
       d3==amap.null_dart_handle || d4==amap.null_dart_handle ||
       d5==amap.null_dart_handle || d6==amap.null_dart_handle ) return false;

  if (!is_face_combinatorial_polygon(amap, d1, 4) ||
      !is_face_combinatorial_polygon(amap, d2, 4) ||
      !is_face_combinatorial_polygon(amap, d3, 4) ||
      !is_face_combinatorial_polygon(amap, d4, 4) ||
      !is_face_combinatorial_polygon(amap, d5, 4) ||
      !is_face_combinatorial_polygon(amap, d6, 4) ) return false;

  // TODO do better with marks.
  if ( belong_to_same_cell<Map,2,1>(amap, d1, d2) ||
       belong_to_same_cell<Map,2,1>(amap, d1, d3) ||
       belong_to_same_cell<Map,2,1>(amap, d1, d4) ||
       belong_to_same_cell<Map,2,1>(amap, d1, d5) ||
       belong_to_same_cell<Map,2,1>(amap, d1, d6) ||
       belong_to_same_cell<Map,2,1>(amap, d2, d3) ||
       belong_to_same_cell<Map,2,1>(amap, d2, d4) ||
       belong_to_same_cell<Map,2,1>(amap, d2, d5) ||
       belong_to_same_cell<Map,2,1>(amap, d2, d6) ||
       belong_to_same_cell<Map,2,1>(amap, d3, d4) ||
       belong_to_same_cell<Map,2,1>(amap, d3, d5) ||
       belong_to_same_cell<Map,2,1>(amap, d3, d6) ||
       belong_to_same_cell<Map,2,1>(amap, d4, d5) ||
       belong_to_same_cell<Map,2,1>(amap, d4, d6) ||
       belong_to_same_cell<Map,2,1>(amap, d5, d6) )
    return false;

  if ( amap.beta(d1,2)    !=amap.beta(d4,1,1) ||
       amap.beta(d1,1,2)  !=amap.beta(d6,0)   ||
       amap.beta(d3,1,2)  !=amap.beta(d6,1)   ||
       amap.beta(d3,0,2)  !=amap.beta(d5,1,1) ||
       amap.beta(d6,1,1,2)!=amap.beta(d2,1)   ||
       amap.beta(d5,0,2)  !=amap.beta(d4,0)   ||
       amap.beta(d5,1,2)  !=amap.beta(d2,0) ) return false;

  return true;
}

  /** Create a new combinatorial hexahedron.
   * @param amap the used combinatorial map.
   * @return a new dart. 
   */
  template < class Map >
  typename Map::Dart_handle make_combinatorial_hexahedron(Map& amap)
  {
    typename Map::Dart_handle d1 = make_combinatorial_polygon(amap,4);
    typename Map::Dart_handle d2 = make_combinatorial_polygon(amap,4);
    typename Map::Dart_handle d3 = make_combinatorial_polygon(amap,4);
    typename Map::Dart_handle d4 = make_combinatorial_polygon(amap,4);
    typename Map::Dart_handle d5 = make_combinatorial_polygon(amap,4);
    typename Map::Dart_handle d6 = make_combinatorial_polygon(amap,4);

    return make_combinatorial_hexahedron(amap, d1, d2, d3, d4, d5, d6);
  }
  
} // namespace CGAL

#endif // CGAL_COMBINATORIAL_MAP_CONSTRUCTORS_H //
// EOF //
