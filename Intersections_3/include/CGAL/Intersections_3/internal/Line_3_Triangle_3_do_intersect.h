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
// Author(s)     : Philippe Guigue

#ifndef CGAL_INTERNAL_INTERSECTIONS_LINE_3_TRIANGLE_3_DO_INTERSECT_H
#define CGAL_INTERNAL_INTERSECTIONS_LINE_3_TRIANGLE_3_DO_INTERSECT_H

#include <CGAL/enum.h>
#include <CGAL/kernel_assertions.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
bool do_intersect(const typename K::Triangle_3& t,
                  const typename K::Line_3& l,
                  const K& k)
{
  CGAL_kernel_precondition(!k.is_degenerate_3_object()(t));
  CGAL_kernel_precondition(!k.is_degenerate_3_object()(l));

  typedef typename K::Point_3 Point_3;

  typename K::Construct_point_on_3 point_on = k.construct_point_on_3_object();
  typename K::Construct_vertex_3 vertex_on = k.construct_vertex_3_object();
  typename K::Orientation_3 orientation = k.orientation_3_object();
  typename K::Coplanar_orientation_3 coplanar_orientation = k.coplanar_orientation_3_object();

  const Point_3& a = vertex_on(t,0);
  const Point_3& b = vertex_on(t,1);
  const Point_3& c = vertex_on(t,2);
  const Point_3& p = point_on(l,0);
  const Point_3& q = point_on(l,1);

  if((orientation(a,b,c,p) != COPLANAR) || (orientation(a,b,c,q) != COPLANAR))
  {
    const Orientation pqab = orientation(p,q,a,b);
    const Orientation pqbc = orientation(p,q,b,c);

    switch(pqab)
    {
      case POSITIVE: return (pqbc != NEGATIVE) && (orientation(p,q,c,a) != NEGATIVE);
      case NEGATIVE: return (pqbc != POSITIVE) && (orientation(p,q,c,a) != POSITIVE);
      case COPLANAR:
        switch(pqbc)
        {
          case POSITIVE: return (orientation(p,q,c,a) != NEGATIVE);
          case NEGATIVE: return (orientation(p,q,c,a) != POSITIVE);
          case COPLANAR: return true;
          default: // should not happen.
            CGAL_kernel_assertion(false);
            return false;
        }
      default: // should not happen.
        CGAL_kernel_assertion(false);
        return false;
    }
  }

  // Coplanar case
  const Orientation pqa = coplanar_orientation(p,q,a);
  return (coplanar_orientation(p,q,b) != pqa) || (coplanar_orientation(p,q,c) != pqa);
}

template <class K>
inline
bool do_intersect(const typename K::Line_3& l,
                  const typename K::Triangle_3& t,
                  const K& k)
{
  return do_intersect(t, l, k);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_LINE_3_TRIANGLE_3_DO_INTERSECT_H
