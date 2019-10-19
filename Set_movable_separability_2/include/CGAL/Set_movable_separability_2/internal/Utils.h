// Copyright (c) 2005,2006,2007,2008,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Shahar    <shasha94@gmail.com>
//            Efi Fogel <efif@gmail.com>

#ifndef CGAL_SET_MOVABLE_SEPARABILITY_2_INTERNAL_UTILS_H
#define CGAL_SET_MOVABLE_SEPARABILITY_2_INTERNAL_UTILS_H

#include <CGAL/license/Set_movable_separability_2.h>


#include <CGAL/enum.h>
#include <CGAL/Polygon_2.h>

namespace CGAL {
namespace Set_movable_separability_2 {
namespace internal {

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
  return (orientation == CGAL::CLOCKWISE) ?
    std::make_pair(backward, forward) : std::make_pair(forward, backward);
}

template <typename Kernel>
bool is_any_edge_colinear(const CGAL::Polygon_2<Kernel>& pgn, Kernel& kernel)
{
  typedef typename Kernel::Point_2                      Point_2;
  typedef typename CGAL::Polygon_2<Kernel>              Polygon_2;
  typedef typename Polygon_2::Vertex_const_iterator     Vertex_const_iterator;
  auto collinear = kernel.collinear_2_object();
  Vertex_const_iterator vci = pgn.vertices_begin();
  Point_2 firstVar = *(vci++);
  Point_2 secondVar = *(vci++);
  Point_2 thirdVar = *(vci++);
  for (; vci != pgn.vertices_end(); ++vci) {
    firstVar = secondVar;
    secondVar = thirdVar;
    thirdVar = *vci;
    if (collinear(firstVar, secondVar, thirdVar)) return true;
  }
  vci = pgn.vertices_begin();
  firstVar = secondVar;
  secondVar = thirdVar;
  thirdVar = *(vci++);
  if(collinear(firstVar, secondVar, thirdVar)) return true;

  firstVar = secondVar;
  secondVar = thirdVar;
  thirdVar = *(vci++);
  if (collinear(firstVar, secondVar, thirdVar)) return true;

  return false;
}

} // namespace internal
} // namespace Set_movable_separability_2
} // namespace CGAL

#endif
