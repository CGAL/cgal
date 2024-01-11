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

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_PLANE_3_PLANE_3_PLANE_3_INTERSECTION_H
#define CGAL_INTERNAL_INTERSECTIONS_3_PLANE_3_PLANE_3_PLANE_3_INTERSECTION_H

#include <CGAL/Intersection_traits_3.h>
#include <CGAL/Intersections_3/internal/Plane_3_Plane_3_intersection.h>
#include <CGAL/Intersections_3/internal/Line_3_Plane_3_intersection.h>

#include <optional>
#include <variant>

namespace CGAL {
namespace Intersections {
namespace internal {

//  triple plane intersection
template <class K>
std::optional<typename K::Point_3>
intersection_point(const typename K::Plane_3& plane1,
                   const typename K::Plane_3& plane2,
                   const typename K::Plane_3& plane3,
                   const K&)
{
  typedef typename K::FT FT;

  const FT& m00 = plane1.a();
  const FT& m01 = plane1.b();
  const FT& m02 = plane1.c();
  const FT& b0  = - plane1.d();
  const FT& m10 = plane2.a();
  const FT& m11 = plane2.b();
  const FT& m12 = plane2.c();
  const FT& b1  = - plane2.d();
  const FT& m20 = plane3.a();
  const FT& m21 = plane3.b();
  const FT& m22 = plane3.c();
  const FT& b2  = - plane3.d();

  // Minors common to two determinants
  const FT minor_0 = m00*m11 - m10*m01;
  const FT minor_1 = m00*m21 - m20*m01;
  const FT minor_2 = m10*m21 - m20*m11;

  const FT den = minor_0*m22 - minor_1*m12 + minor_2*m02; // determinant of M

  if(is_zero(den)){
    return std::nullopt;
  }

  const FT num3 = minor_0*b2 - minor_1*b1 + minor_2*b0;  // determinant of M with M[x:2] swapped with [b0,b1,b2]

  // Minors common to two determinants
  const FT minor_3 = b0*m12 - b1*m02;
  const FT minor_4 = b0*m22 - b2*m02;
  const FT minor_5 = b1*m22 - b2*m12;

  // num1 has opposite signs because b0 and M[:1] have been swapped
  const FT num1 = - minor_3*m21 + minor_4*m11 - minor_5*m01;  // determinant of M with M[x:0] swapped with [b0,b1,b2]
  const FT num2 = minor_3*m20 - minor_4*m10 + minor_5*m00;  // determinant of M with M[x:1] swapped with [b0,b1,b2]

  return std::make_optional(typename K::Point_3(num1/den, num2/den, num3/den));
}

template <class K>
std::optional<std::variant<typename K::Point_3, typename K::Line_3, typename K::Plane_3> >
intersection(const typename K::Plane_3& plane1,
             const typename K::Plane_3& plane2,
             const typename K::Plane_3& plane3,
             const K& k)
{
  typedef typename std::optional<std::variant<typename K::Point_3,
                                                  typename K::Line_3,
                                                  typename K::Plane_3> > result_type;

  typedef typename K::Line_3       Line_3;
  typedef typename K::Plane_3      Plane_3;

  auto res = intersection_point(plane1,plane2,plane3, k);
  if (res)
    return result_type(*res);

  // Intersection between plane1 and plane2 can either be
  // a line, a plane, or empty.
  typename Intersection_traits<K, Plane_3, Plane_3>::result_type
      o12 = internal::intersection(plane1, plane2, k);

  if(o12)
  {
    if(const Line_3* l = intersect_get<Line_3>(o12))
    {
      if (internal::do_intersect(*l, plane3, k))
        return result_type(*l);
    }
    else
    {
      CGAL_assertion(intersect_get<Plane_3>(o12) != nullptr);
      // either line or plane
      typename Intersection_traits<K, Plane_3, Plane_3>::result_type
          v = internal::intersection(plane3, plane1, k);
      if(v)
      {
        if( intersect_get<Plane_3>(v)!=nullptr)
          return result_type(plane1);
        else if(const Line_3* l = intersect_get<Line_3>(v))
          return result_type(*l);
      }
    }
  }

  return result_type();
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_PLANE_3_PLANE_3_PLANE_3_INTERSECTION_H
