// Copyright (c) 2024 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_ISOSURFACING_3_INTERNAL_IMPLICIT_SHAPES_HELPER_H
#define CGAL_ISOSURFACING_3_INTERNAL_IMPLICIT_SHAPES_HELPER_H

#include <CGAL/number_utils.h>
#include <CGAL/Kernel/global_functions_3.h>

namespace CGAL {
namespace Isosurfacing {
namespace Shapes {

// Shapes are defined at isovalue 0

// c is the center
// r the radius
template<typename K>
typename K::FT
sphere(const typename K::Point_3& c,
       const typename K::FT r,
       const typename K::Point_3& q)
{
  return CGAL::approximate_sqrt(CGAL::squared_distance(c, q)) - r;
}

template <typename K>
typename K::FT
box(const typename K::Point_3& b,
    const typename K::Point_3& t,
    const typename K::Point_3& q)
{
  typename K::Point_3 c = CGAL::midpoint(b, t);
  typename K::Iso_cuboid_3 ic(b, t);
  bool inside = ic.has_on_bounded_side(q);
  typename K::FT d = 0;
  if(inside)
  {
    d = (std::min)({CGAL::abs(q.x() - b.x()), CGAL::abs(q.x() - t.x()),
                    CGAL::abs(q.y() - b.y()), CGAL::abs(q.y() - t.y()),
                    CGAL::abs(q.z() - b.z()), CGAL::abs(q.z() - t.z())});
  }
  else
  {
    for(int i=0; i<3; ++i)
      d += (CGAL::abs(q[i] - c[i]) > (c[i] - b[i]) ? CGAL::square(q[i] - c[i]) : 0);
    d = CGAL::approximate_sqrt(d);
  }

  return inside ? - d : d;
}

// template <typename K>
// typename K::FT
// disk(const typename K::Point_3& c,
//      const typename K::Vector_3& n,
//      const typename K::FT r,
//      const typename K::Point_3& q)
// {
//   typename K::Plane_3 pl(c, n);
//   typename K::Point_3 pq = pl.projection(q);

//   typename K::FT sq_dpl = CGAL::squared_distance(q, pl);

//   if(CGAL::squared_distance(pq, c) < CGAL::square(r))
//     return CGAL::approximate_sqrt(sq_dpl);
//   else
//     return CGAL::approximate_sqrt(CGAL::square(CGAL::approximate_sqrt(CGAL::squared_distance(pq, c)) - r) + sq_dpl);
// }



// p is the center of the base disk
// q is the center of the top disk
template<typename K>
typename K::FT
infinite_cylinder(const typename K::Point_3& b,
                  const typename K::Vector_3& n,
                  const typename K::FT r,
                  const typename K::Point_3& q)
{
  typename K::Plane_3 pl(b, n);
  typename K::Point_3 pq = pl.projection(q);
  return CGAL::approximate_sqrt(CGAL::squared_distance(pq, b)) - r;
}


// c is the center of the torus
// n is the normal of the plane containing all centers of the tube
// r is the small radius
// R is the large radius
template<typename K>
typename K::FT
torus(const typename K::Point_3& c,
      const typename K::Vector_3& n,
      const typename K::FT r,
      const typename K::FT R,
      const typename K::Point_3& q)
{
  typename K::Vector_3 w (c, q);
  typename K::Plane_3 pl(c, n);
  typename K::Point_3 pq = pl.projection(q);
  typename K::FT d = CGAL::approximate_sqrt(CGAL::squared_distance(pq, c)) - R;
  typename K::FT h = CGAL::abs(CGAL::scalar_product(w, n));

  return CGAL::approximate_sqrt(CGAL::square(d) + CGAL::square(h)) - r;
}

template<typename K>
typename K::FT
torus_ridge(const typename K::Point_3& c,
            const typename K::Vector_3& n,
            const typename K::FT r,
            const typename K::FT R,
            const typename K::Point_3& q)
{
  typename K::Vector_3 w (c, q);
  typename K::Plane_3 pl(c, n);
  typename K::Point_3 pq = pl.projection(q);
  typename K::FT d = CGAL::abs(CGAL::approximate_sqrt(CGAL::squared_distance(pq, c)) - R) - r;
  return d + CGAL::squared_distance(q, pl);
}

template<typename K>
typename K::FT
inverted_torus(const typename K::Point_3& c,
               const typename K::Vector_3& n,
               const typename K::FT r,
               const typename K::FT R,
               const typename K::Point_3& q)
{
  typename K::Vector_3 w (c, q);
  typename K::Plane_3 pl(c, n);
  typename K::Point_3 pq = pl.projection(q);
  typename K::FT d = CGAL::abs(CGAL::approximate_sqrt(CGAL::squared_distance(pq, c)) - R) - r;
  return d - CGAL::squared_distance(q, pl);
}

/////////////////////////////////////////////////////////////////

template <typename K, typename S1, typename S2>
typename K::FT
shape_union(const S1& s1, const S2& s2, const typename K::Point_3& q)
{
  return std::min(s1(q), s2(q));
}

template <typename K, typename S1, typename S2>
typename K::FT
shape_difference(const S1& s1, const S2& s2, const typename K::Point_3& q)
{
  return std::max(s1(q), -s2(q));
}

template <typename K, typename S1, typename S2>
typename K::FT
shape_intersection(const S1& s1, const S2& s2, const typename K::Point_3& q)
{
  return std::max(s1(q), s2(q));
}

template <typename K, typename S1, typename S2>
typename K::FT
shape_symmetric_difference(const S1& s1, const S2& s2, const typename K::Point_3& q)
{
  return std::max(-std::min(s1(q), s2(q)), std::max(s1(q), s2(q)));
}


} // namespace Shapes
} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_INTERNAL_IMPLICIT_SHAPES_HELPER_H
