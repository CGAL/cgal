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
// Author(s)     : Sébastien Loriot

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_RAY_3_SEGMENT_3_DO_INTERSECT_H
#define CGAL_INTERNAL_INTERSECTIONS_3_RAY_3_SEGMENT_3_DO_INTERSECT_H

#include <CGAL/Intersections_3/internal/Line_3_Segment_3_do_intersect.h>
#include <CGAL/Intersections_3/internal/Point_3_Ray_3_do_intersect.h>

#include <CGAL/enum.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
inline
bool
do_intersect(const typename K::Segment_3& s,
             const typename K::Ray_3& r,
             const K& k)
{
  CGAL_precondition(!s.is_degenerate() && !r.is_degenerate());

  if(!do_intersect(s, r.supporting_line(), k))
    return false;

  typename K::Coplanar_orientation_3 pred = k.coplanar_orientation_3_object();

  CGAL::Orientation p0p1s = pred(s.point(0), s.point(1), r.source());
  CGAL::Orientation stp0 = pred(r.source(), r.second_point(), s.point(0));

  if(p0p1s == COLLINEAR) //s belongs to the supporting line of p0p1
  {
    if(stp0 == COLLINEAR )//st and p0p1 have the same supporting line
      return Ray_3_has_on_collinear_Point_3(r, s.point(0), k) ||
             Ray_3_has_on_collinear_Point_3(r, s.point(1), k);
    else
      return true;
  }

  if(stp0 == COLLINEAR)
    return Ray_3_has_on_collinear_Point_3(r, s.point(0), k);

  return (p0p1s != stp0);
}

template <class K>
inline
bool
do_intersect(const typename K::Ray_3& r,
             const typename K::Segment_3& s,
             const K& k)
{
  return do_intersect(s, r, k);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_RAY_3_SEGMENT_3_DO_INTERSECT_H
