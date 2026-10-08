// Copyright (c) 2019
// GeometryFactory.  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Maxime Gimeno

#ifndef CGAL_INTERSECTIONS_2_CIRCLE_2_SEGMENT_2_H
#define CGAL_INTERSECTIONS_2_CIRCLE_2_SEGMENT_2_H

#include <CGAL/Distance_2/internal/squared_distance_utils_2.h>
#include <CGAL/Intersection_traits_2.h>
#include <CGAL/tags.h>

#include <CGAL/Circle_2.h>
#include <CGAL/Segment_2.h>

namespace CGAL {
namespace Intersections {
namespace internal {
//


template <class K>
typename K::Boolean
do_intersect(const typename K::Circle_2& circ,
             const typename K::Segment_2& s,
             const K& k)
{
  typedef typename K::Vector_2 Vector_2;
  typedef typename K::Point_2 Point_2;
  typedef typename K::FT FT;

  typename K::Compute_scalar_product_2 sp = k.compute_scalar_product_2_object();
  typename K::Compute_squared_distance_2 sq_dist = k.compute_squared_distance_2_object();
  typename K::Compute_squared_length_2 sq_length = k.compute_squared_length_2_object();
  typename K::Compute_squared_radius_2 squared_radius = k.compute_squared_radius_2_object();
  typename K::Construct_center_2 center = k.construct_center_2_object();
  typename K::Construct_line_2 line = k.construct_line_2_object();
  typename K::Construct_source_2 source = k.construct_source_2_object();
  typename K::Construct_target_2 target = k.construct_target_2_object();
  typename K::Construct_vector_2 vector = k.construct_vector_2_object();

  const Point_2 c = center(circ);
  const FT sqrad = squared_radius(circ);
  const Vector_2 src_dir = vector(source(s), c);
  const Vector_2 tgt_dir = vector(target(s), c);
  const Vector_2 seg_dir = vector(source(s), target(s));

  const FT src_sp = sp(src_dir, seg_dir);
  const FT tgt_sp = sp(tgt_dir, seg_dir);

  // Is segment pointing away from circle?
  if (!is_positive(src_sp)) {
    FT lsrc = sq_length(src_dir);
    if (lsrc > sqrad)
      return false;

    FT ltgt = sq_length(tgt_dir);
    return ltgt >= sqrad;
  }

  // Is segment pointing towards circle?
  if (!is_negative(tgt_sp)) {
    FT ltgt = sq_length(tgt_dir);
    if (ltgt > sqrad)
      return false;

    FT lsrc = sq_length(src_dir);
    return lsrc >= sqrad;
  }

  FT lsrc = sq_length(src_dir);
  FT ltgt = sq_length(tgt_dir);

  // Fully contained in circle
  if (lsrc < sqrad && ltgt < sqrad)
      return false;

  return sq_dist(c, line(s)) <= sqrad;
}

template <class K>
typename K::Boolean
do_intersect(const typename K::Segment_2& s,
             const typename K::Circle_2& c,
             const K& k)
{
  return do_intersect(c, s, k);
}

} // namespace internal
} // namespace Intersections

CGAL_DO_INTERSECT_FUNCTION(Circle_2, Segment_2, 2)

} // namespace CGAL

#endif // CGAL_INTERSECTIONS_2_CIRCLE_2_SEGMENT_2_H
