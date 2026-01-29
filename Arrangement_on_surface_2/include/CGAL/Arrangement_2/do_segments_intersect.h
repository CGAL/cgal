// Copyright (c) 202Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_AOS_2_DO_SEGMENTS_INTERSECT_H
#define CGAL_AOS_2_DO_SEGMENTS_INTERSECT_H

#include <CGAL/license/Arrangement_on_surface_2.h>

/*! \file Intersection detection functions.
 */

namespace CGAL {
namespace Aos_2 {
namespace internal {

// Specialized do_intersect with several tests skipped because at
// this point, we already know the order of points
template <typename Point_2, typename Traits>
bool do_closed_segment_intersect(const Point_2& l1, const Point_2& r1, const Point_2& l2, const Point_2& r2,
                                 const Traits& traits) {
  auto cmpare_xy = traits.compare_xy_2_object();
  auto compute_orientation = traits.orientation_2_object();

  auto res1 = make_certain(cmpare_xy(l1, l2));
  if (res1 == EQUAL) return true;
  auto res2 = make_certain(cmpare_xy(r1, r2));
  if (res2 == EQUAL) return true;

  if (res1 == SMALLER) {
    // (2, 3) l1 < l2
    auto orient1 = compute_orientation(l1, r1, l2);
    if (res2 == SMALLER) {
      // (2) r1 < r2
      if (orient1 == COLLINEAR) return true;
      auto orient2 = compute_orientation(l2, r2, r1);
      return ((orient2 == COLLINEAR) || (orient2 == orient1));
    }
    // (3) r2 < r1
    if (orient1 == COLLINEAR) return true;
    auto orient2 = compute_orientation(l1, r1, r2);
    return ((orient2 == COLLINEAR) || (orient2 != orient1));
  }

  // (5, 6) l2 < r1
  auto orient1 = compute_orientation(l2, r2, l1);
  if (res2 == SMALLER) {
    // (6) r1 < r2
    if (orient1 == COLLINEAR) return true;
    auto orient2 = compute_orientation(l2, r2, r1);
    return ((orient2 == COLLINEAR) || (orient2 != orient1));
  }
  // (5) r2 < r1
  if (orient1 == COLLINEAR) return true;
  auto orient2 = compute_orientation(l1, r1, r2);
  return ((orient2 == COLLINEAR) || (orient2 == orient1));
}

// Specialized do_intersect with several tests skipped because at
// this point, we already know the order of points
template <typename Point_2, typename Traits>
bool do_open_segment_intersect(const Point_2& l1, const Point_2& r1, const Point_2& l2, const Point_2& r2,
                               const Traits& traits) {
  auto cmpare_xy = traits.compare_xy_2_object();
  auto compute_orientation = traits.orientation_2_object();

  auto res1 = make_certain(cmpare_xy(l1, l2));
  if (res1 == EQUAL) return (make_certain(compute_orientation(l1, r1, r2)) == COLLINEAR);
  auto res2 = make_certain(cmpare_xy(r1, r2));
  if (res2 == EQUAL) return (make_certain(compute_orientation(r2, l2, l1)) == COLLINEAR);

  if (res1 == SMALLER) {
    // (2, 3) l1 < l2
    auto orient1 = compute_orientation(l1, r1, l2);
    if (res2 == SMALLER) {
      // (2) r1 < r2
      if (orient1 == COLLINEAR) return true;
      auto orient2 = compute_orientation(l2, r2, r1);
      return ((orient2 == COLLINEAR) || (orient2 == orient1));
    }
    // (3) r2 < r1
    if (orient1 == COLLINEAR) return true;
    auto orient2 = compute_orientation(l1, r1, r2);
    return ((orient2 == COLLINEAR) || (orient2 != orient1));
  }

  // (5, 6) l2 < l1
  auto orient1 = compute_orientation(l2, r2, l1);
  if (res2 == SMALLER) {
    // (6) r1 < r2
    if (orient1 == COLLINEAR) return true;
    auto orient2 = compute_orientation(l2, r2, r1);
    return ((orient2 == COLLINEAR) || (orient2 != orient1));
  }
  // (5) r2 < r1
  if (orient1 == COLLINEAR) return true;
  auto orient2 = compute_orientation(l1, r1, r2);
  return ((orient2 == COLLINEAR) || (orient2 == orient1));
}

/*! determines whether two segments intersect.
 */
template <typename Segment, typename Traits>
do_segment_intersect(const Segment& seg1, const Segment& seg2, bool closed, const Traits& traits) {
  auto ctr_min_vertex = traits.construct_min_vertex_2_object();
  auto ctr_max_vertex = traits.construct_max_vertex_2_object();
  auto cmpare_xy = traits.compare_xy_2_object();

  const auto& l1 = ctr_min_vertex(seg1);
  const auto& r1 = ctr_max_vertex(seg1);
  const auto& l2 = ctr_min_vertex(seg2);
  const auto& r2 = ctr_max_vertex(seg2);

  // There are 6 cases (not including degenerate cases):
  // 1. l1, r1, l2, r2
  // 2. l1, l2, r1, r2
  // 3. l1, l2, r2, r1
  // 4. l2, r2, l1, r1
  // 5. l2, l1, r2, r1
  // 6. l2, l1, r1, r2
  // Handle cases (1) and (4) first (relatively trivial):
  switch (make_certain(cmpare_xy(r1, l2))) {
   case SMALLER: return false;
   case EQUAL: return closed;
   default: break; // LERGER
  }
  switch (make_certain(cmpare_xy(r2, l1))) {
   case SMALLER: return false;
   case EQUAL: return closed;
   default: break; // LERGER
  }

  return (closed) ? do_closed_segment_intersect(l1, r1, l2, r2, traits) : do_open_segment_intersect(l1, r1, l2, r2, traits);
}

}
}
}
#endif
