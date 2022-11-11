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
// Author(s)     : Mariette Yvinec

#ifndef CGAL_INTERNAL_INTERSECTIONS_POINT_3_TRIANGLE_3_DO_INTERSECT_H
#define CGAL_INTERNAL_INTERSECTIONS_POINT_3_TRIANGLE_3_DO_INTERSECT_H

#include <CGAL/enum.h>
#include <CGAL/kernel_assertions.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
bool do_intersect(const typename K::Triangle_3& t,
                  const typename K::Point_3& p,
                  const K& k)
{
  CGAL_kernel_precondition(!k.is_degenerate_3_object()(t));

  typedef typename K::Point_3 Point_3;

  typename K::Construct_vertex_3 vertex_on = k.construct_vertex_3_object();
  typename K::Orientation_3 orientation = k.orientation_3_object();
  typename K::Coplanar_orientation_3 coplanar_orientation = k.coplanar_orientation_3_object();

  const Point_3& a = vertex_on(t,0);
  const Point_3& b = vertex_on(t,1);
  const Point_3& c = vertex_on(t,2);

  if(orientation(a,b,c,p) != COPLANAR)
    return false;

  const CGAL::Orientation abp = coplanar_orientation(a,b,p);
  const CGAL::Orientation bcp = coplanar_orientation(b,c,p);

  switch(abp)
  {
  case POSITIVE:
    return (bcp != NEGATIVE) && (coplanar_orientation(c,a,p) != NEGATIVE);
  case NEGATIVE:
    return (bcp != POSITIVE) && (coplanar_orientation(c,a,p) != POSITIVE);
  case COLLINEAR:
    switch (bcp)
    {
    case POSITIVE:
      return (coplanar_orientation(c,a,p) != NEGATIVE);
    case NEGATIVE:
      return (coplanar_orientation(c,a,p) != POSITIVE);
    case COLLINEAR: return true;
    default: // should not happen.
      CGAL_kernel_assertion(false);
      return false;
    }
  default: // should not happen.
    CGAL_kernel_assertion(false);
    return false;
  }
}

template <class K>
bool do_intersect(const typename K::Point_3& p,
                  const typename K::Triangle_3& t,
                  const K& k)
{
  return do_intersect(t, p, k);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_POINT_3_TRIANGLE_3_DO_INTERSECT_H
