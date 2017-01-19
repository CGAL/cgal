// Copyright (c) 2015  Tel-Aviv University (Israel).
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
// Author(s): Sebastian Morr    <sebastian@morr.cc>

#ifndef CGAL_MINKOWSKI_SUM_HOLE_FILTER_2_H
#define CGAL_MINKOWSKI_SUM_HOLE_FILTER_2_H

#include <CGAL/license/Minkowski_sum_2.h>


#include <CGAL/basic.h>
#include <vector>

namespace CGAL {

/*! \class
 * This class applies filter to a polygon with holes,
 * by removing all of its holes that cannot possibly contribute
 * to the Minkowski sum boundary.
 */
template <typename Kernel_, typename Container_>
class Hole_filter_2
{
private:
  typedef Kernel_               Kernel;
  typedef Container_            Container;

  typedef CGAL::Polygon_2<Kernel, Container>            Polygon_2;
  typedef CGAL::Polygon_with_holes_2<Kernel, Container> Polygon_with_holes_2;
  typedef typename Polygon_with_holes_2::Hole_iterator  Hole_iterator;
  typedef std::vector<Hole_iterator>                    Hole_iterator_vector;

public:
  /*! Filter out holes of a polygon with holes.
   * \param[in] pgn1 The polygon with holes to filter.
   * \param[in] pgn2 The reference polygon with holes.
   * \param[out] filtered_pgn1 the filterd polygon.
   */
  void operator()(const Polygon_with_holes_2& pgn1,
                  const Polygon_2& pgn2,
                  Polygon_with_holes_2& filtered_pgn1) const
  {
    filtered_pgn1 = pgn1;

    Hole_iterator_vector to_erase;
    Bbox_2 boundary_bbox = pgn2.bbox();

    Hole_iterator it = filtered_pgn1.holes_begin();
    while (it != filtered_pgn1.holes_end()) {
      Bbox_2 hole_bbox = (*it).bbox();

      if ((hole_bbox.ymax()-hole_bbox.ymin() <
           boundary_bbox.ymax()-boundary_bbox.ymin()) ||
          (hole_bbox.xmax()-hole_bbox.xmin() <
           boundary_bbox.xmax()-boundary_bbox.xmin()))
      {
        to_erase.push_back(it);
      }
      ++it;
    }

    typename Hole_iterator_vector::iterator it2 = to_erase.begin();
    while (it2 != to_erase.end()) filtered_pgn1.erase_hole(*it2++);
  }

  /*! Filter out holes of a polygon with holes.
   * \param[in] pgn1 The polygon with holes to filter.
   * \param[in] pgn2 The reference polygon polygon with holes.
   * \param[out] filtered_pgn1 the filterd polygon.
   */
  void operator()(const Polygon_with_holes_2& pgn1,
                  const Polygon_with_holes_2& pgn2,
                  Polygon_with_holes_2& filtered_pgn1) const
  { operator()(pgn1, pgn2.outer_boundary(), filtered_pgn1); }
};

} // namespace CGAL

#endif
