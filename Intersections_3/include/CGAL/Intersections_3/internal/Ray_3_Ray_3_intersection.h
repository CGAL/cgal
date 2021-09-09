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

#ifndef CGAL_RAY_3_RAY_3_INTERSECTION_H
#define CGAL_RAY_3_RAY_3_INTERSECTION_H

#include <CGAL/Intersections_3/internal/Line_3_Ray_3_intersection.h>

#include <CGAL/enum.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
typename Intersection_traits<K, typename K::Ray_3, typename K::Ray_3>::result_type
intersection(const typename K::Ray_3& r1,
             const typename K::Ray_3& r2,
             const K& k)
{
  CGAL_precondition(!r1.is_degenerate() && !r2.is_degenerate());

  typename Intersection_traits<K, typename K::Line_3, typename K::Ray_3>::result_type
      v = internal::intersection(r1.supporting_line(), r2, k);

  if(v)
  {
    if(const typename K::Point_3* p = intersect_get<typename K::Point_3>(v))
    {
      if(Ray_3_has_on_collinear_Point_3(r1,*p,k))
        return intersection_return<typename K::Intersect_3, typename K::Ray_3, typename K::Ray_3>(*p);
    }
    else if(const typename K::Ray_3* r = intersect_get<typename K::Ray_3>(v))
    {
      bool r1_has_s2 = Ray_3_has_on_collinear_Point_3(r1,r2.source(),k);
      bool r2_has_s1 = Ray_3_has_on_collinear_Point_3(r2,r1.source(),k);

      if(r1_has_s2)
      {
        if(r2_has_s1)
        {
          if(k.equal_3_object()(r1.source(),r2.source()))
          {
            return intersection_return<typename K::Intersect_3, typename K::Ray_3, typename K::Ray_3>(r1.source());
          }
          else
          {
            return intersection_return<typename K::Intersect_3, typename K::Ray_3, typename K::Ray_3>(
                     k.construct_segment_3_object()(r1.source(), r2.source()));
          }
        }
        else
        {
          return intersection_return<typename K::Intersect_3, typename K::Ray_3, typename K::Ray_3>(*r);
        }
      }
      else
      {
        if(r2_has_s1)
          return intersection_return<typename K::Intersect_3, typename K::Ray_3, typename K::Ray_3>(r1);
      }
    }
  }

  return intersection_return<typename K::Intersect_3, typename K::Ray_3, typename K::Ray_3>();
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_RAY_3_RAY_3_INTERSECTION_H
