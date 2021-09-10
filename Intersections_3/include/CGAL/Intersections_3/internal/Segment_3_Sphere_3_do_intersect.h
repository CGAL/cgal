// Copyright (c) 2003  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_INTERNAL_INTERSECTIONS_SEGMENT_3_SPHERE_3_DO_INTERSECT_H
#define CGAL_INTERNAL_INTERSECTIONS_SEGMENT_3_SPHERE_3_DO_INTERSECT_H

#include <CGAL/Rational_traits.h>
#include <CGAL/Distance_3/Point_3_Segment_3.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
inline
typename K::Boolean
do_intersect(const typename K::Sphere_3& sp,
             const typename K::Segment_3& seg,
             const K& k)
{
  typedef typename K::RT RT;
  typedef typename K::Bounded_side Bounded_side;

  typename K::Bounded_side_3 bounded_side = k.bounded_side_3_object();
  typename K::Compute_squared_radius_3 sq_radius = k.compute_squared_radius_3_object();
  typename K::Construct_center_3 center = k.construct_center_3_object();
  typename K::Construct_source_3 source = k.construct_source_3_object();
  typename K::Construct_target_3 target = k.construct_target_3_object();

  const Bounded_side source_side = bounded_side(sp, source(seg));
  const Bounded_side target_side = bounded_side(sp, target(seg));

  if(source_side != target_side)
    return true;
  else if(source_side == ON_BOUNDED_SIDE) // else if ==> both extremities are on the same side
    return false;
  else if(source_side == ON_BOUNDARY)
    return true;

  CGAL_kernel_assertion(source_side == ON_UNBOUNDED_SIDE && target_side == ON_UNBOUNDED_SIDE);

  // Both extremities are out, but the segment can still be intersecting the sphere
  RT num, den;
  CGAL::internal::squared_distance_RT(center(sp), seg, num, den, k);

  return !(compare_quotients<RT>(num, den,
                                 Rational_traits<typename K::FT>().numerator(sq_radius(sp)),
                                 Rational_traits<typename K::FT>().denominator(sq_radius(sp))) == LARGER);
}

template <class K>
inline
typename K::Boolean
do_intersect(const typename K::Segment_3& seg,
             const typename K::Sphere_3& sp,
             const K& k)
{
  return do_intersect(sp, seg, k);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_SEGMENT_3_SPHERE_3_DO_INTERSECT_H
