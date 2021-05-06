// Copyright (c) 2018 GeometryFactory Sarl
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_TRIANGLE_3_SPHERE_3_DO_INTERSECT_H
#define CGAL_TRIANGLE_3_SPHERE_3_DO_INTERSECT_H

#include <CGAL/squared_distance_3_2.h>
#include <CGAL/Intersection_traits_3.h>
#include <CGAL/Rational_traits.h>
namespace CGAL {

template <class K>
class Triangle_3;

template <class K>
class Sphere_3;

template <class K>
class Line_3;

namespace Intersections {

namespace internal {

template <class K>
inline
typename K::Boolean
do_intersect(const typename K::Sphere_3 &sp,
             const typename K::Triangle_3 &tr,
             const K & k)
{
  typedef typename K::RT RT;
  RT num, den;

  CGAL::internal::squared_distance_RT(sp.center(), tr, num, den, k);
  return ! (compare_quotients<RT>(num, den,
                                  Rational_traits<typename K::FT>().numerator(sp.squared_radius()),
                                  Rational_traits<typename K::FT>().denominator(sp.squared_radius())) == LARGER);
}

template <class K>
inline
typename K::Boolean
do_intersect(const typename K::Triangle_3 &tr,
             const typename K::Sphere_3 &sp,
             const K & k)
{
  return do_intersect(sp, tr, k);
}


template <class K>
inline
typename K::Boolean
do_intersect(const typename K::Line_3 &lin,
             const typename K::Sphere_3 &sp,
             const K & k)
{
  typedef typename K::RT RT;
  RT num, den;

  CGAL::internal::squared_distance_RT(sp.center(), lin, num, den, k);
  return ! (compare_quotients<RT>(num, den,
                                  Rational_traits<typename K::FT>().numerator(sp.squared_radius()),
                                  Rational_traits<typename K::FT>().denominator(sp.squared_radius())) == LARGER);

}

template <class K>
inline
typename K::Boolean
do_intersect(const typename K::Sphere_3 &sp,
             const typename K::Line_3 &lin,
             const K & k)
{
  return do_intersect(lin,sp,k);
}



template <class K>
inline
typename K::Boolean
do_intersect(const typename K::Sphere_3 &sp,
             const typename K::Ray_3 &ray,
             const K & k)
{
  typedef typename K::RT RT;
  RT num, den;

  CGAL::internal::squared_distance_RT(sp.center(), ray, num, den, k);
  return ! (compare_quotients<RT>(num, den,
                                  Rational_traits<typename K::FT>().numerator(sp.squared_radius()),
                                  Rational_traits<typename K::FT>().denominator(sp.squared_radius())) == LARGER);
}


template <class K>
inline
typename K::Boolean
do_intersect(const typename K::Ray_3 &ray,
             const typename K::Sphere_3 &sp,
             const K & k)
{
  return do_intersect(sp,ray,k);
}

template <class K>
inline
typename K::Boolean
do_intersect(const typename K::Sphere_3 &sp,
             const typename K::Segment_3 &seg,
             const K & k)
{
  typedef typename K::RT RT;
  RT num, den;

  CGAL::internal::squared_distance_RT(sp.center(), seg, num, den, k);
  return ! (compare_quotients<RT>(num, den,
                                  Rational_traits<typename K::FT>().numerator(sp.squared_radius()),
                                  Rational_traits<typename K::FT>().denominator(sp.squared_radius())) == LARGER);
}


template <class K>
inline
typename K::Boolean
do_intersect(const typename K::Segment_3 &seg,
             const typename K::Sphere_3 &sp,
             const K & k)
{
  return do_intersect(sp,seg,k);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_TRIANGLE_3_SPHERE_3_DO_INTERSECT_H
