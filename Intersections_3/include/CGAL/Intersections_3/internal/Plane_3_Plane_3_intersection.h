// Copyright (c) 2008 INRIA(France).
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

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_PLANE_3_PLANE_3_INTERSECTION_H
#define CGAL_INTERNAL_INTERSECTIONS_3_PLANE_3_PLANE_3_INTERSECTION_H

#include <CGAL/Intersection_traits_3.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
typename Intersection_traits<K, typename K::Plane_3, typename K::Plane_3>::result_type
intersection(const typename K::Plane_3& plane1,
             const typename K::Plane_3& plane2,
             const K&)
{
  typedef typename K::RT RT;
  typedef typename K::Point_3 Point_3;
  typedef typename K::Direction_3 Direction_3;
  typedef typename K::Line_3 Line_3;

  const RT& a = plane1.a();
  const RT& b = plane1.b();
  const RT& c = plane1.c();
  const RT& d = plane1.d();
  const RT& p = plane2.a();
  const RT& q = plane2.b();
  const RT& r = plane2.c();
  const RT& s = plane2.d();

  RT det = a*q-p*b;
  if(det != 0)
  {
    Point_3 is_pt = Point_3(b*s-d*q, p*d-a*s, 0, det);
    Direction_3 is_dir = Direction_3(b*r-c*q, p*c-a*r, det);
    return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Plane_3>(Line_3(is_pt, is_dir));
  }

  det = a*r-p*c;
  if(det != 0)
  {
    Point_3 is_pt = Point_3(c*s-d*r, 0, p*d-a*s, det);
    Direction_3 is_dir = Direction_3(c*q-b*r, det, p*b-a*q);
    return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Plane_3>(Line_3(is_pt, is_dir));
  }

  det = b*r-c*q;
  if(det != 0)
  {
    Point_3 is_pt = Point_3(0, c*s-d*r, d*q-b*s, det);
    Direction_3 is_dir = Direction_3(det, c*p-a*r, a*q-b*p);
    return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Plane_3>(Line_3(is_pt, is_dir));
  }

  // degenerate case
  if(a!=0 || p!=0)
  {
    if(a*s == p*d)
      return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Plane_3>(plane1);
    else
      return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Plane_3>();
  }

  if(b!=0 || q!=0)
  {
    if(b*s == q*d)
      return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Plane_3>(plane1);
    else
      return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Plane_3>();
  }

  if(c!=0 || r!=0)
  {
    if(c*s == r*d)
      return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Plane_3>(plane1);
    else
      return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Plane_3>();
  }

  return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Plane_3>(plane1);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_PLANE_3_PLANE_3_INTERSECTION_H
