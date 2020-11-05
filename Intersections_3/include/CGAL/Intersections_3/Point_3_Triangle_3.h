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

#ifndef CGAL_INTERSECTIONS_3_POINT_3_TRIANGLE_3_H
#define CGAL_INTERSECTIONS_3_POINT_3_TRIANGLE_3_H

#include <CGAL/Triangle_3.h>
#include <CGAL/Point_3.h>

#include <CGAL/enum.h>
#include <CGAL/kernel_assertions.h>
#include <CGAL/Intersection_traits_3.h>

namespace CGAL {

namespace Intersections {

namespace internal {

template <class K>
bool do_intersect(const typename K::Triangle_3 &t,
                  const typename K::Point_3    &p,
                  const K & k )
{

  CGAL_kernel_precondition( ! k.is_degenerate_3_object()(t));

  typedef typename K::Point_3 Point_3;

  typename K::Construct_vertex_3 vertex_on =
    k.construct_vertex_3_object();

  typename K::Orientation_3 orientation =
    k.orientation_3_object();

  typename K::Coplanar_orientation_3 coplanar_orientation =
    k.coplanar_orientation_3_object();



  const Point_3 & a = vertex_on(t,0);
  const Point_3 & b = vertex_on(t,1);
  const Point_3 & c = vertex_on(t,2);


  if (orientation(a,b,c,p) != COPLANAR)
    return false;


  const Orientation abp = coplanar_orientation(a,b,p);
  const Orientation bcp = coplanar_orientation(b,c,p);


  switch ( abp ) {
  case POSITIVE: return  bcp != NEGATIVE
                   &&   coplanar_orientation(c,a,p) != NEGATIVE ;
  case NEGATIVE: return  bcp != POSITIVE
                   &&   coplanar_orientation(c,a,p) != POSITIVE ;
  case COLLINEAR:
    switch ( bcp ) {
    case POSITIVE: return  coplanar_orientation(c,a,p) != NEGATIVE ;
    case NEGATIVE: return  coplanar_orientation(c,a,p) != POSITIVE ;
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
bool do_intersect(const typename K::Point_3    &p,
                  const typename K::Triangle_3 &t,
                  const K & k )
{
  return do_intersect(t, p, k);
}

template <class K>
inline
typename CGAL::Intersection_traits
<K, typename K::Point_3, typename K::Triangle_3>::result_type
intersection(const typename K::Point_3 &pt,
             const typename K::Triangle_3 &tr,
             const K& k)
{
  if (do_intersect(pt,tr, k)) {
    return intersection_return<typename K::Intersect_3, typename K::Point_3, typename K::Triangle_3>(pt);
  }
  return intersection_return<typename K::Intersect_3, typename K::Point_3, typename K::Triangle_3>();
}

template <class K>
inline
typename CGAL::Intersection_traits
<K, typename K::Triangle_3, typename K::Point_3>::result_type
intersection( const typename K::Triangle_3 &tr,
              const typename K::Point_3 &pt,
              const K& k)
{
  return internal::intersection(pt, tr, k);
}
} // namespace internal
} // namespace Intersections

CGAL_DO_INTERSECT_FUNCTION(Triangle_3, Point_3, 3)
CGAL_INTERSECTION_FUNCTION(Triangle_3, Point_3, 3)

} //namespace CGAL

#endif // CGAL_INTERSECTIONS_3_POINT_3_TRIANGLE_3_H
