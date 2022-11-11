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

#ifndef CGAL_INTERNAL_INTERSECTIONS_RAY_3_RAY_3_DO_INTERSECT_H
#define CGAL_INTERNAL_INTERSECTIONS_RAY_3_RAY_3_DO_INTERSECT_H

#include <CGAL/Intersections_3/internal/Line_3_Ray_3_do_intersect.h>
#include <CGAL/Intersections_3/internal/Point_3_Ray_3_do_intersect.h>

#include <CGAL/enum.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
inline
bool
do_intersect(const typename K::Ray_3& r1,
             const typename K::Ray_3& r2,
             const K& k)
{
  CGAL_precondition(!r1.is_degenerate() && !r2.is_degenerate());

  if(!do_intersect(r1, r2.supporting_line(), k))
    return false;

  typename K::Coplanar_orientation_3 pred = k.coplanar_orientation_3_object();

  CGAL::Orientation p0p1s = pred(r1.source(), r1.second_point(), r2.source());
  CGAL::Orientation stp0 = pred(r2.source(), r2.second_point(), r1.source());

  if(p0p1s == COLLINEAR)
  {
    if(stp0 == COLLINEAR)
      return Ray_3_has_on_collinear_Point_3(r2, r1.source(), k) ||
             Ray_3_has_on_collinear_Point_3(r1, r2.source(), k);
    else
      return true;
  }

  if(stp0 == COLLINEAR)
    return Ray_3_has_on_collinear_Point_3(r2, r1.source(), k);

  return (p0p1s != stp0);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_RAY_3_RAY_3_DO_INTERSECT_H
