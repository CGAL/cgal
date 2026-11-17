// Copyright (c) 2026 GeometryFactory Sarl
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sébastien Loriot
//

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_BALL_3_BALL_3_DO_INTERSECT_H
#define CGAL_INTERNAL_INTERSECTIONS_3_BALL_3_BALL_3_DO_INTERSECT_H

#include <CGAL/number_utils.h>

namespace CGAL {
namespace Intersections {
namespace internal {

// Two balls intersect iff the distance between their centers is not larger
// than the sum of their radii, i.e. d2 <= sr1 + sr2 + 2*sqrt(sr1*sr2),
// which is evaluated without square root.
template <class K>
inline
typename K::Boolean
do_intersect(const typename K::Ball_3& b1,
             const typename K::Ball_3& b2,
             const K& k)
{
  typedef typename K::FT FT;

  typename K::Construct_center_3 center = k.construct_center_3_object();
  typename K::Compute_squared_radius_3 sqr = k.compute_squared_radius_3_object();
  typename K::Compute_squared_distance_3 sq_dist = k.compute_squared_distance_3_object();

  const FT sr1 = sqr(b1);
  const FT sr2 = sqr(b2);
  const FT s = sq_dist(center(b1), center(b2)) - sr1 - sr2;

  // we want to know if d(c1,c2) <= r1+r2

  // r1^2+r2^2 is enough to cover d(c1,c2)^2
  if(s <= FT(0))
    return true;

  // squaring again the whole thing
  return CGAL::square(s) <= FT(4) * sr1 * sr2;
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_BALL_3_BALL_3_DO_INTERSECT_H
