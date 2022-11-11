// Copyright (c) 2018 INRIA Sophia-Antipolis (France)
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

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_LINE_3_RAY_3_DO_INTERSECT_H
#define CGAL_INTERNAL_INTERSECTIONS_3_LINE_3_RAY_3_DO_INTERSECT_H

#include <CGAL/Intersections_3/internal/Point_3_Ray_3_do_intersect.h>

#include <CGAL/enum.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
inline
bool
do_intersect(const typename K::Line_3& l,
             const typename K::Ray_3& r,
             const K& k)
{
  CGAL_precondition(!l.is_degenerate() && !r.is_degenerate());

  if(!do_intersect(l, r.supporting_line()))
    return false;

  typename K::Coplanar_orientation_3 pred = k.coplanar_orientation_3_object();
  CGAL::Orientation p0p1s = pred(l.point(0), l.point(1), r.source());
  if(p0p1s == COLLINEAR)
    return true;

  CGAL::Orientation stp0 = pred(r.source(), r.second_point(), l.point(0));
  if(stp0 == COLLINEAR)
    return Ray_3_has_on_collinear_Point_3(r,l.point(0),k);

  return (p0p1s != stp0);
}

template <class K>
inline
bool
do_intersect(const typename K::Ray_3& r,
             const typename K::Line_3& l,
             const K& k)
{
  return do_intersect(l,r,k);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_LINE_3_RAY_3_DO_INTERSECT_H
