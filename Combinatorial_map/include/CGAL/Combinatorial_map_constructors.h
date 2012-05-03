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
    amap.template basic_link_beta<2>(d1, d2);
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
      amap.template basic_link_beta<1>(prev, cur);
      prev=cur;
    }
  
    amap.template basic_link_beta<1>(prev, start);
    return start;
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
    amap.basic_link_beta(d1, d2, 2);
    amap.basic_link_beta(d3, d2->beta(0), 2);
    amap.basic_link_beta(d1->beta(1), d3->beta(0), 2);
    amap.basic_link_beta(d4, d2->beta(1), 2);
    amap.basic_link_beta(d4->beta(0), d3->beta(1), 2);
    amap.basic_link_beta(d4->beta(1), d1->beta(0), 2);

    return d1;
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
    amap.basic_link_beta(d1,                   d4->beta(1)->beta(1), 2);
    amap.basic_link_beta(d1->beta(1),          d6->beta(0)         , 2);
    amap.basic_link_beta(d1->beta(1)->beta(1), d2                  , 2);
    amap.basic_link_beta(d1->beta(0),          d5                  , 2);

    amap.basic_link_beta(d3,                   d2->beta(1)->beta(1), 2);
    amap.basic_link_beta(d3->beta(1),          d6->beta(1)         , 2);
    amap.basic_link_beta(d3->beta(1)->beta(1), d4                  , 2);
    amap.basic_link_beta(d3->beta(0),          d5->beta(1)->beta(1), 2);

    amap.basic_link_beta(d6,                   d4->beta(1)         , 2);
    amap.basic_link_beta(d6->beta(1)->beta(1), d2->beta(1)         , 2);

    amap.basic_link_beta(d5->beta(0),          d4->beta(0)         , 2);
    amap.basic_link_beta(d5->beta(1),          d2->beta(0)         , 2);

    return d1;
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
