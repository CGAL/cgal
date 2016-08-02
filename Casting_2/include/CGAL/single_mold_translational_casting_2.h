// Copyright (c) 2005,2006,2007,2008,2009,2010,2011 Tel-Aviv University (Israel).
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
// Author(s): Shahar    <shasha94@gmail.com>
//            Efi Fogel <efif@gmail.com>

#ifndef CGAL_SINGLE_MOLD_TRANSLATIONAL_CASTING_2_H
#define CGAL_SINGLE_MOLD_TRANSLATIONAL_CASTING_2_H

#include <CGAL/Polygon_2.h>
#include <CGAL/enum.h>

#include <iostream>
#include <list>
#include "Casting_2/Circle_arrangment.h"
namespace CGAL {

/* Legend:
 * point = Represented as Direction_2. It is the intersection between the
 *   fitting Direction_2 and the unit circle
 *
 * Arc = Represented as A pair of point. clockwise arc between the first
 *   point and the second point. (each of its sides might be open or closed)
 *
 * SegmentOuterCircle  = Arc that represent all the directions that points
 *   out from the polygon if it start from the
 *   fitting segment. This arc is always open half circle.
 */

/*! \fn std::pair<typename Kernel::Direction_2,typename Kernel::Direction_2> get_segment_outer_circle(typename Kernel::Segment_2 seg, CGAL::Orientation orientation)
 * \param[in] seg the polygon segment
 * \param[in] orientation the orientation of the segment (and the polygon).
 *   if CLOCKWISE then the outer half circle is to the left.
 * \return the open outer half-circle of the edge.
 */
  template <typename Kernel>
  inline std::pair<typename Kernel::Direction_2, typename Kernel::Direction_2>
  get_segment_outer_circle(const typename Kernel::Segment_2 seg,
                           const CGAL::Orientation orientation)
  {
    typename Kernel::Direction_2 forward( seg);
    typename Kernel::Direction_2 backward(-forward);
    return (orientation == CGAL::Orientation::CLOCKWISE) ?
      std::make_pair(backward, forward) : std::make_pair(forward, backward);
  }

  template <typename Kernel>
  bool isAnyEdgeColinear(const CGAL::Polygon_2<Kernel>& pgn)
  {
    typedef typename CGAL::Point_2<Kernel>        Point_2;
    typedef typename CGAL::Polygon_2<Kernel>        Polygon_2;
    typedef typename Polygon_2::Vertex_const_iterator     Vertex_const_iterator;
    Vertex_const_iterator vci = pgn.vertices_begin();
    Point_2 firstVar  = *(vci++);
    Point_2 secondVar  = *(vci++);
    Point_2 thirdVar  = *(vci++);
    for(;vci!=pgn.vertices_end();++vci)
    {
      firstVar=secondVar;
      secondVar=thirdVar;
      thirdVar = *vci;
      if(CGAL::collinear(firstVar,secondVar,thirdVar)) return true;
    }
    vci = pgn.vertices_begin();
    firstVar=secondVar;
    secondVar=thirdVar;
    thirdVar = *(vci++);
    if(CGAL::collinear(firstVar,secondVar,thirdVar)) return true;

    firstVar=secondVar;
    secondVar=thirdVar;
    thirdVar = *(vci++);
    if(CGAL::collinear(firstVar,secondVar,thirdVar)) return true;

    return false;
  }


  /*! \fn OutputIterator find_single_mold_translational_casting_2(const CGAL::Polygon_2<Kernel>& pgn, OutputIterator oi)
   * \param[in] pgn the input polygon that we want to check if is castable or not.
   * \param[in,out] oi the output iterator to put the top edges in
   * \return all the possible top edges of the polygon and there pullout direction
   *   (with no rotation)
   */
  template <typename Kernel, typename OutputIterator>
  OutputIterator
  single_mold_translational_casting_2(const CGAL::Polygon_2<Kernel>& pgn,
                                      OutputIterator oi)
  {
    /* Legend
     * point = Represented as  Direction_2. It is the intersection between the
     *   fitting Direction_2 and the unit circle
     *
     * arc = Represented as A pair of point. clockwise arc between the first
     *   point and the second point. (each of its sides might be open or closed)
     */
    CGAL_precondition(pgn.is_simple());
    CGAL_precondition(!isAnyEdgeColinear(pgn));


    auto e_it = pgn.edges_begin();
    size_t edge_index = 0;
    CGAL::Orientation poly_orientation = pgn.orientation();
    auto segment_outer_circle =
      get_segment_outer_circle<Kernel>(*e_it++, poly_orientation);
    INNER_CASTING_2::Circle_arrangment<Kernel>
      circle_arrangment(segment_outer_circle);

    ++edge_index;
    for (; e_it!= pgn.edges_end(); ++e_it,++edge_index) {
      segment_outer_circle =
        get_segment_outer_circle<Kernel>(*e_it, poly_orientation);
      circle_arrangment.add_segment_outer_circle(segment_outer_circle, edge_index);
      if (circle_arrangment.all_is_covered_twice()) return oi;
    }
    circle_arrangment.get_all_1_edges(oi);
    return oi;
  }


}

#endif
