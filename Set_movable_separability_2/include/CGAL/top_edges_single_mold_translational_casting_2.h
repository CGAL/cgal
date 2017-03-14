// Copyright (c) 2016 Tel-Aviv University (Israel).
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

#ifndef CGAL_TOP_EDGES_SINGLE_MOLD_TRANSLATIONAL_CASTING_2_H
#define CGAL_TOP_EDGES_SINGLE_MOLD_TRANSLATIONAL_CASTING_2_H

#include <iostream>
#include <list>

#include <CGAL/Polygon_2.h>
#include <CGAL/enum.h>

#include "Set_movable_separability_2/Circle_arrangment.h"
#include "Set_movable_separability_2/Utils.h"

namespace CGAL {

  namespace Set_movable_separability_2 {

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

    /*! \fn OutputIterator find_single_mold_translational_casting_2(const CGAL::Polygon_2<Kernel>& pgn, OutputIterator oi)
     * \param[in] pgn the input polygon that we want to check if is castable or not.
     * \param[in,out] oi the output iterator to put the top edges in
     * \param[in] kernel the kernel to use.
     * \return all the possible top edges of the polygon and there pullout direction
     *  a pair of Directions is build this way [firstClockwise,secondClockwise]
     *   (with no rotation)
     */
    template <typename Kernel, typename OutputIterator>
    OutputIterator
    top_edges_single_mold_translational_casting_2
    (const CGAL::Polygon_2<Kernel>& pgn, OutputIterator oi, Kernel& kernel)
    {
      /* Legend
       * point = Represented as  Direction_2. It is the intersection between the
       *   fitting Direction_2 and the unit circle
       *
       * arc = Represented as A pair of point. clockwise arc between the first
       *   point and the second point. (each of its sides might be open or closed)
       */
      CGAL_precondition(pgn.is_simple());
      CGAL_precondition(!internal::is_any_edge_colinear(pgn));

      auto e_it = pgn.edges_begin();
      CGAL::Orientation poly_orientation = pgn.orientation();
      auto segment_outer_circle =
	  internal::get_segment_outer_circle<Kernel>(*e_it++, poly_orientation);
      internal::Circle_arrangment<Kernel> circle_arrangment(kernel,
							    segment_outer_circle,pgn.edges_begin());

      for (; e_it != pgn.edges_end(); ++e_it) {
	  segment_outer_circle =
	      internal::get_segment_outer_circle<Kernel>(*e_it, poly_orientation);
	  circle_arrangment.add_segment_outer_circle(segment_outer_circle, e_it);
	  if (circle_arrangment.all_is_covered_twice()) return oi;
      }
      circle_arrangment.get_all_1_edges(oi);
      return oi;
    }

    /*! \fn OutputIterator find_single_mold_translational_casting_2(const CGAL::Polygon_2<Kernel>& pgn, OutputIterator oi)
     * \param[in] pgn the input polygon that we want to check if is castable or not.
     * \param[in,out] oi the output iterator to put the top edges in
     * \return all the possible top edges of the polygon and there pullout direction
     *  a pair of Directions is build this way [firstClockwise,secondClockwise]
     *   (with no rotation)
     */
    template <typename Kernel, typename OutputIterator>
    OutputIterator
    top_edges_single_mold_translational_casting_2
    (const CGAL::Polygon_2<Kernel>& pgn, OutputIterator oi)
    {
      Kernel kernel;
      return top_edges_single_mold_translational_casting_2(pgn, oi, kernel);
    }

  } // end of namespace Set_movable_separability_2
} // end of namespace CGAL

#endif
