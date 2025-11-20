// Copyright (c) 2025 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef PRISM_AND_PYRAMID_CREATION_H
#define PRISM_AND_PYRAMID_CREATION_H

/** Create a combinatorial prism from 2 triangles and 3 squares.
   * @param amap the used combinatorial map.
   * @param d1 a dart onto a first triangle.
   * @param d2 a dart onto a first square.
   * @param d3 a dart onto a second square.
   * @param d4 a dart onto a thirth square.
   * @param d5 a dart onto a second triangle.
   * @return a new dart.
   */
  template < class Map >
  typename Map::Dart_handle
  make_combinatorial_prism(Map& amap,
                           typename Map::Dart_handle d1,
                           typename Map::Dart_handle d2,
                           typename Map::Dart_handle d3,
                           typename Map::Dart_handle d4,
                           typename Map::Dart_handle d5)
  {
    // 2-link for first triangle
    amap.basic_link_beta_for_involution(d1, d2, 2);
    amap.basic_link_beta_for_involution(amap.beta(d1, 1), d3, 2);
    amap.basic_link_beta_for_involution(amap.beta(d1, 0), d4, 2);

    // 2-link for quandrangles between them
    amap.basic_link_beta_for_involution(amap.beta(d2, 0), amap.beta(d3, 1), 2);
    amap.basic_link_beta_for_involution(amap.beta(d2, 1), amap.beta(d4, 0), 2);
    amap.basic_link_beta_for_involution(amap.beta(d3, 0), amap.beta(d4, 1), 2);

    // 2-link for second triangle
    amap.basic_link_beta_for_involution(amap.beta(d2, 1, 1), d5, 2);
    amap.basic_link_beta_for_involution(amap.beta(d3, 1, 1), amap.beta(d5, 0), 2);
    amap.basic_link_beta_for_involution(amap.beta(d4, 1, 1), amap.beta(d5, 1), 2);

    return d1;
  }

  /** Create a new combinatorial prism.
   * @param amap the used combinatorial map.
   * @return a new dart.
   */
  template < class Map >
  typename Map::Dart_handle make_combinatorial_prism(Map& amap)
  {
    typename Map::Dart_handle d1 = make_combinatorial_polygon(amap,3);
    typename Map::Dart_handle d2 = make_combinatorial_polygon(amap,4);
    typename Map::Dart_handle d3 = make_combinatorial_polygon(amap,4);
    typename Map::Dart_handle d4 = make_combinatorial_polygon(amap,4);
    typename Map::Dart_handle d5 = make_combinatorial_polygon(amap,3);

    return make_combinatorial_prism(amap, d1, d2, d3, d4, d5);
  }

/** Test if a volume is a combinatorial prism.
 * @param amap the used combinatorial map.
 * @param adart an intial dart
 * @return true iff the volume containing adart is a combinatorial prism.
 */
template < class Map >
bool is_volume_combinatorial_prism(const Map& amap,
                                   typename Map::Dart_const_handle d1)
{
  typename Map::Dart_const_handle d2 = amap.beta(d1, 2);
  typename Map::Dart_const_handle d3 = amap.beta(d1, 1, 2);
  typename Map::Dart_const_handle d4 = amap.beta(d1, 0, 2);
  typename Map::Dart_const_handle d5 = amap.beta(d2, 1, 1, 2);

  if ( d1==amap.null_dart_handle || d2==amap.null_dart_handle ||
       d3==amap.null_dart_handle || d4==amap.null_dart_handle ||
       d5==amap.null_dart_handle ) return false;

  if (!amap.is_face_combinatorial_polygon(d1, 3) ||
      !amap.is_face_combinatorial_polygon(d2, 4) ||
      !amap.is_face_combinatorial_polygon(d3, 4) ||
      !amap.is_face_combinatorial_polygon(d4, 4) ||
      !amap.is_face_combinatorial_polygon(d5, 3)) return false;

  // TODO do better with marks.
  if ( amap.template belong_to_same_cell<2,1>(d1, d2) ||
       amap.template belong_to_same_cell<2,1>(d1, d3) ||
       amap.template belong_to_same_cell<2,1>(d1, d4) ||
       amap.template belong_to_same_cell<2,1>(d1, d5) ||
       amap.template belong_to_same_cell<2,1>(d2, d3) ||
       amap.template belong_to_same_cell<2,1>(d2, d4) ||
       amap.template belong_to_same_cell<2,1>(d2, d5) ||
       amap.template belong_to_same_cell<2,1>(d3, d4) ||
       amap.template belong_to_same_cell<2,1>(d3, d5) ||
       amap.template belong_to_same_cell<2,1>(d4, d5))
    return false;

  if (amap.beta(d2,0,2)  !=amap.beta(d3,1) ||
      amap.beta(d2,1,2)  !=amap.beta(d4,0) ||
      amap.beta(d3,0,2)  !=amap.beta(d4,1) ||
      amap.beta(d3,1,1,2)!=amap.beta(d5,0) ||
      amap.beta(d4,1,1,2)!=amap.beta(d5,1)) return false;

  return true;
}

/** Create a combinatorial pyramid from 1 square and 4 triangles.
   * @param amap the used combinatorial map.
   * @param d1 a dart onto the square.
   * @param d2 a dart onto a first triangle.
   * @param d3 a dart onto a second triangle.
   * @param d4 a dart onto a thirth triangle.
   * @param d5 a dart onto a fourth triangle.
   * @return a new dart.
   */
  template < class Map >
  typename Map::Dart_handle
  make_combinatorial_pyramid(Map& amap,
                             typename Map::Dart_handle d1,
                             typename Map::Dart_handle d2,
                             typename Map::Dart_handle d3,
                             typename Map::Dart_handle d4,
                             typename Map::Dart_handle d5)
  {
    // 2-link for the square
    amap.basic_link_beta_for_involution(d1, d2, 2);
    amap.basic_link_beta_for_involution(amap.beta(d1, 1), d5, 2);
    amap.basic_link_beta_for_involution(amap.beta(d1, 1, 1), d4, 2);
    amap.basic_link_beta_for_involution(amap.beta(d1, 0), d3, 2);

    // 2-link for first triangle
    amap.basic_link_beta_for_involution(amap.beta(d2, 1), amap.beta(d3, 0), 2);
    amap.basic_link_beta_for_involution(amap.beta(d2, 0), amap.beta(d5, 1), 2);

    // 2-link for triangles between them
    amap.basic_link_beta_for_involution(amap.beta(d5, 0), amap.beta(d4, 1), 2);
    amap.basic_link_beta_for_involution(amap.beta(d4, 0), amap.beta(d3, 1), 2);

    return d1;
  }

  /** Create a new combinatorial pyramid.
   * @param amap the used combinatorial map.
   * @return a new dart.
   */
  template < class Map >
  typename Map::Dart_handle make_combinatorial_pyramid(Map& amap)
  {
    typename Map::Dart_handle d1 = make_combinatorial_polygon(amap,4);
    typename Map::Dart_handle d2 = make_combinatorial_polygon(amap,3);
    typename Map::Dart_handle d3 = make_combinatorial_polygon(amap,3);
    typename Map::Dart_handle d4 = make_combinatorial_polygon(amap,3);
    typename Map::Dart_handle d5 = make_combinatorial_polygon(amap,3);

    return make_combinatorial_pyramid(amap, d1, d2, d3, d4, d5);
  }

/** Test if a volume is a combinatorial pyramid.
 * @param amap the used combinatorial map.
 * @param adart an intial dart
 * @return true iff the volume containing adart is a combinatorial pyramid.
 */
template < class Map >
bool is_volume_combinatorial_pyramid(const Map& amap,
                                     typename Map::Dart_const_handle d1)
{
  typename Map::Dart_const_handle d2 = amap.beta(d1, 2);
  typename Map::Dart_const_handle d3 = amap.beta(d1, 0, 2);
  typename Map::Dart_const_handle d4 = amap.beta(d1, 1, 1, 2);
  typename Map::Dart_const_handle d5 = amap.beta(d1, 1, 2);

  if ( d1==amap.null_dart_handle || d2==amap.null_dart_handle ||
       d3==amap.null_dart_handle || d4==amap.null_dart_handle ||
       d5==amap.null_dart_handle ) return false;

  if (!amap.is_face_combinatorial_polygon(d1, 4) ||
      !amap.is_face_combinatorial_polygon(d2, 3) ||
      !amap.is_face_combinatorial_polygon(d3, 3) ||
      !amap.is_face_combinatorial_polygon(d4, 3) ||
      !amap.is_face_combinatorial_polygon(d5, 3)) return false;

  // TODO do better with marks.
  if ( amap.template belong_to_same_cell<2,1>(d1, d2) ||
       amap.template belong_to_same_cell<2,1>(d1, d3) ||
       amap.template belong_to_same_cell<2,1>(d1, d4) ||
       amap.template belong_to_same_cell<2,1>(d1, d5) ||
       amap.template belong_to_same_cell<2,1>(d2, d3) ||
       amap.template belong_to_same_cell<2,1>(d2, d4) ||
       amap.template belong_to_same_cell<2,1>(d2, d5) ||
       amap.template belong_to_same_cell<2,1>(d3, d4) ||
       amap.template belong_to_same_cell<2,1>(d3, d5) ||
       amap.template belong_to_same_cell<2,1>(d4, d5))
    return false;

  if (amap.beta(d2,1,2)!=amap.beta(d3,0) ||
      amap.beta(d2,0,2)!=amap.beta(d5,1) ||
      amap.beta(d5,0,2)!=amap.beta(d4,1) ||
      amap.beta(d4,0,2)!=amap.beta(d3,1)) return false;

  return true;
}

/** Create an prism given 6 Vertex_attribute_handle.
     *    (6 vertices, 9 edges and 5 facets)
     * \verbatim
     *      3---4
     *      |\ /|
     *      0-5-1
     *       \|/
     *        2
     * \endverbatim
     * @param h0 the first vertex handle.
     * @param h1 the second vertex handle.
     * @param h2 the third vertex handle.
     * @param h3 the fourth vertex handle.
     * @param h4 the fifth vertex handle.
     * @param h5 the sixth vertex handle.
     * @return the dart of the new prism incident to h0 and to
     *         the facet (h0,h1,h2).
     */
     template < class Map >
     typename Map::Dart_handle make_prism(Map& amap,
                                          typename Map::Vertex_attribute_handle h0,
                                          typename Map::Vertex_attribute_handle h1,
                                          typename Map::Vertex_attribute_handle h2,
                                          typename Map::Vertex_attribute_handle h3,
                                          typename Map::Vertex_attribute_handle h4,
                                          typename Map::Vertex_attribute_handle h5)
    {
      typename Map::Dart_handle d1 = amap.make_triangle(h0, h1, h2);
      typename Map::Dart_handle d2 = amap.make_quadrangle(h1, h0, h3, h4);
      typename Map::Dart_handle d3 = amap.make_quadrangle(h2, h1, h4, h5);
      typename Map::Dart_handle d4 = amap.make_quadrangle(h0, h2, h5, h3);
      typename Map::Dart_handle d5 = amap.make_triangle(h4, h3, h5);

      return make_combinatorial_prism(amap, d1, d2, d3, d4, d5);
    }

    /** Create an prism given 6 points.
     * \verbatim
     *      3---4
     *      |\ /|
     *      0-5-1
     *       \|/
     *        2
     * \endverbatim
     * @param p0 the first point.
     * @param p1 the second point.
     * @param p2 the third point.
     * @param p3 the fourth point.
     * @param p4 the fifth point.
     * @param p5 the sixth point.
     * @return the dart of the new prism incident to p0 and to
     *         the facet (p0,p1,p2).
     */
     template < class Map >
     typename Map::Dart_handle make_prism(Map& amap,
                            const typename Map::Point& p0,
                            const typename Map::Point& p1,
                            const typename Map::Point& p2,
                            const typename Map::Point& p3,
                            const typename Map::Point& p4,
                            const typename Map::Point& p5)
    {
      return make_prism(amap,
                        amap.create_vertex_attribute(p0),
                        amap.create_vertex_attribute(p1),
                        amap.create_vertex_attribute(p2),
                        amap.create_vertex_attribute(p3),
                        amap.create_vertex_attribute(p4),
                        amap.create_vertex_attribute(p5));
    }

    /** Create an pyramid given 5 Vertex_attribute_handle.
     *    (5 vertices, 8 edges and 5 facets)
     * \verbatim
     *       4
     *      /|\
     *     0-|-1
     *     | | |
     *     3---2
     * \endverbatim
     * @param h0 the first vertex handle.
     * @param h1 the second vertex handle.
     * @param h2 the third vertex handle.
     * @param h3 the fourth vertex handle.
     * @param h4 the fifth vertex handle.
     * @return the dart of the new pyramid incident to h0 and to
     *         the facet (h0,h1,h2,h3).
     */
     template < class Map >
    typename Map::Dart_handle make_pyramid(Map& amap,
                             typename Map::Vertex_attribute_handle h0,
                             typename Map::Vertex_attribute_handle h1,
                             typename Map::Vertex_attribute_handle h2,
                             typename Map::Vertex_attribute_handle h3,
                             typename Map::Vertex_attribute_handle h4)
    {
      typename Map::Dart_handle d1 = amap.make_quadrangle(h0, h1, h2, h3);
      typename Map::Dart_handle d2 = amap.make_triangle(h1, h0, h4);
      typename Map::Dart_handle d3 = amap.make_triangle(h0, h3, h4);
      typename Map::Dart_handle d4 = amap.make_triangle(h3, h2, h4);
      typename Map::Dart_handle d5 = amap.make_triangle(h2, h1, h4);

      return make_combinatorial_pyramid(amap, d1, d2, d3, d4, d5);
    }

    /** Create an pyramid given 5 points.
     * \verbatim
     *       4
     *      /|\
     *     0-|-1
     *     | | |
     *     3---2
     * \endverbatim
     * @param p0 the first point.
     * @param p1 the second point.
     * @param p2 the third point.
     * @param p3 the fourth point.
     * @param p4 the fifth point.
     * @return the dart of the new pyramid incident to p0 and to
     *         the facet (p0,p1,p2,p3).
     */
    template < class Map >
    typename Map::Dart_handle make_pyramid(Map& amap,
                            const typename Map::Point& p0,
                            const typename Map::Point& p1,
                            const typename Map::Point& p2,
                            const typename Map::Point& p3,
                            const typename Map::Point& p4)
    {
      return make_pyramid(amap,
                          amap.create_vertex_attribute(p0),
                          amap.create_vertex_attribute(p1),
                          amap.create_vertex_attribute(p2),
                          amap.create_vertex_attribute(p3),
                          amap.create_vertex_attribute(p4));
    }

#endif // PRISM_AND_PYRAMID_CREATION_H
///////////////////////////////////////////////////////////////////////////////
