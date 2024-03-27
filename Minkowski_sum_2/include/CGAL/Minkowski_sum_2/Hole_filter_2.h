// Copyright (c) 2015  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
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
class Hole_filter_2 {
private:
  using Kernel = Kernel_;
  using Container = Container_;
  using Polygon_2 = CGAL::Polygon_2<Kernel, Container>;
  using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<Kernel, Container>;

public:
  /*! Filter out holes of a polygon with holes.
   * \param[in] pgn1 The polygon with holes to filter.
   * \param[in] pgn2 The reference polygon with holes.
   * \param[out] filtered_pgn1 the filtered polygon.
   */
  void operator()(const Polygon_with_holes_2& pgn1,
                  const Polygon_2& pgn2,
                  Polygon_with_holes_2& filtered_pgn1) const {
    Bbox_2 boundary_bbox = pgn2.bbox();
    filtered_pgn1 = Polygon_with_holes_2(pgn1.outer_boundary());
    for (const auto& h : pgn1.holes()) {
      Bbox_2 hole_bbox = h.bbox();
      if ((hole_bbox.ymax()-hole_bbox.ymin() >=
           boundary_bbox.ymax()-boundary_bbox.ymin()) &&
          (hole_bbox.xmax()-hole_bbox.xmin() >=
           boundary_bbox.xmax()-boundary_bbox.xmin()))
        filtered_pgn1.add_hole(h);
    }
  }

  /*! Filter out holes of a polygon with holes.
   * \param[in] pgn1 The polygon with holes to filter.
   * \param[in] pgn2 The reference polygon polygon with holes.
   * \param[out] filtered_pgn1 the filtered polygon.
   */
  void operator()(const Polygon_with_holes_2& pgn1,
                  const Polygon_with_holes_2& pgn2,
                  Polygon_with_holes_2& filtered_pgn1) const
  { operator()(pgn1, pgn2.outer_boundary(), filtered_pgn1); }
};

} // namespace CGAL

#endif
