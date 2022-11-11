// Copyright (c) 2008 INRIA Sophia-Antipolis (France), ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Geert-Jan Giezeman,
//                 Sylvain Pion

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_LINE_3_PLANE_3_DO_INTERSECT_H
#define CGAL_INTERNAL_INTERSECTIONS_3_LINE_3_PLANE_3_DO_INTERSECT_H

#include <CGAL/wmult.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
bool
do_intersect(const typename K::Plane_3& plane,
             const typename K::Line_3& line,
             const K&)
{
  typedef typename K::Point_3 Point_3;
  typedef typename K::Direction_3 Direction_3;
  typedef typename K::RT RT;

  const Point_3& line_pt = line.point();
  const Direction_3& line_dir = line.direction();

  RT den = plane.a()*line_dir.dx() + plane.b()*line_dir.dy() + plane.c()*line_dir.dz();
  if(den != 0)
    return true;

  RT num = plane.a()*line_pt.hx() + plane.b()*line_pt.hy()
           + plane.c()*line_pt.hz() + wmult_hw((K*)0, plane.d(), line_pt);

  if(num == RT(0)) // all line
    return true;
  else // no intersection
    return false;
}

template <class K>
inline
bool
do_intersect(const typename K::Line_3& line,
             const typename K::Plane_3& plane,
             const K& k)
{
  return do_intersect(plane, line, k);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_LINE_3_PLANE_3_DO_INTERSECT_H
