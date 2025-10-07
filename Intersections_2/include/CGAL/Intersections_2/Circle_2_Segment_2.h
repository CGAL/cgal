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
             const K& k,
             const Homogeneous_tag&)
{
  typedef typename K::Vector_2 Vector_2;
  typedef typename K::Point_2 Point_2;
  typedef typename K::FT FT;

  typename K::Construct_vector_2 vector = k.construct_vector_2_object();
  typename K::Compute_squared_length_2 sq_length = k.compute_squared_length_2_object();
  typename K::Compute_squared_distance_2 sq_dist = k.compute_squared_distance_2_object();

  const Point_2 c = circ.center();
  const FT sqrad = circ.squared_radius();
  const Vector_2 diffsrc = vector(s.source(), c);
  const Vector_2 difftgt = vector(s.target(), c);

  FT lsrc = CGAL::internal::wdot(diffsrc, diffsrc, k);
  FT lsrcw = square(diffsrc.hw());
  FT ltgt = CGAL::internal::wdot(difftgt, difftgt, k);
  FT ltgtw = square(difftgt.hw());

  if (lsrc < sqrad * lsrcw) {
    if (ltgt < sqrad * ltgtw)
      return false;
    else
      return true;
  }
  else if (ltgt <= sqrad * ltgtw)
    return true;

  return sq_dist(c, s.supporting_line()) <= sqrad;
}

template <class K>
typename K::Boolean
do_intersect(const typename K::Circle_2& circ,
             const typename K::Segment_2& s,
             const K& k,
             const Cartesian_tag&)
{
  typedef typename K::Vector_2 Vector_2;
  typedef typename K::Point_2 Point_2;
  typedef typename K::FT FT;

  typename K::Construct_vector_2 vector = k.construct_vector_2_object();
  typename K::Compute_squared_length_2 sq_length = k.compute_squared_length_2_object();
  typename K::Compute_squared_distance_2 sq_dist = k.compute_squared_distance_2_object();

  const Point_2 c = circ.center();
  const FT sqrad = circ.squared_radius();
  const Vector_2 diffsrc = vector(s.source(), c);
  const Vector_2 difftgt = vector(s.target(), c);

  FT lsrc = CGAL::internal::wdot(diffsrc, diffsrc, k);
  FT ltgt = CGAL::internal::wdot(difftgt, difftgt, k);

  if (lsrc < sqrad) {
    if (ltgt < sqrad)
      return false;
    else
      return true;
  }
  else if (ltgt <= sqrad)
    return true;

  return sq_dist(c, s.supporting_line()) <= sqrad;
}

template <class K>
typename K::Boolean
do_intersect(const typename K::Circle_2& c,
             const typename K::Segment_2& s,
             const K& k)
{
  typedef typename K::Kernel_tag Tag;
  Tag tag;
  return do_intersect(c, s, k, tag);
}

template <class K>
typename K::Boolean
do_intersect(const typename K::Segment_2& s,
             const typename K::Circle_2& c,
             const K& k)
{
  typedef typename K::Kernel_tag Tag;
  Tag tag;
  return do_intersect(c, s, k, tag);
}

} // namespace internal
} // namespace Intersections

CGAL_DO_INTERSECT_FUNCTION(Circle_2, Segment_2, 2)

} // namespace CGAL

#endif // CGAL_INTERSECTIONS_2_CIRCLE_2_SEGMENT_2_H
