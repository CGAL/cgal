// Copyright (c) 2009   INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Manuel Caroli <Manuel.Caroli@sophia.inria.fr>

#ifndef CGAL_PERIODIC_3_TRIANGULATION_REMOVE_TRAITS_3_H
#define CGAL_PERIODIC_3_TRIANGULATION_REMOVE_TRAITS_3_H

#include <CGAL/license/Periodic_3_triangulation_3.h>


#include <CGAL/basic.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Periodic_3_offset_3.h>

namespace CGAL {

// Triangulation_3 has calls to Construct_point_3 to handle weighted and bare points.
// The default inherited Construct_point_3 inherited by Periodic_3_triangulation_remove_traits_3
// must be overwritten by a construction Construct_point_3 that offers:
// - pair<K::Point_3, offset> --> pair<K::Point_3, offset> (identity)
// - pair<K::Weighted_point_3, offset> --> pair<K::Point_3, offset> (will be added when needed)
template < class Traits_ >
class Construct_point_3_remove_traits
{
  typedef Traits_                        Traits;
  typedef typename Traits::Point_3       Point_3; // actually a pair <K::Point_3, offset>

public:
  const Point_3& operator()(const Point_3& p) { return p; }
};

template < class Traits_, class Functor_ >
class Point_offset_adaptor {
  typedef Traits_ Traits;
  typedef Functor_ Functor;

  typedef typename Traits::Point_3       Point;

public:
  typedef typename Functor::result_type result_type;

 Point_offset_adaptor(const Functor & functor) : _functor(functor) {}

  result_type operator()(const Point& p0, const Point& p1) const {
    return _functor(p0.first, p1.first,
        p0.second, p1.second);
  }
  result_type operator()(const Point& p0, const Point& p1,
      const Point& p2) const {
    return _functor(p0.first, p1.first, p2.first,
        p0.second, p1.second, p2.second);
  }
  result_type operator()(const Point& p0, const Point& p1,
      const Point& p2, const Point& p3) const {
    return _functor(p0.first, p1.first, p2.first, p3.first,
        p0.second, p1.second, p2.second, p3.second);
  }
  result_type operator()(const Point& p0, const Point& p1,
      const Point& p2, const Point& p3, const Point& p4) const {
    return _functor(p0.first, p1.first, p2.first, p3.first, p4.first,
        p0.second, p1.second, p2.second, p3.second, p4.second);
  }

private:
  Functor _functor;
};

template < class P3DTTraits_, class Off_ = typename CGAL::Periodic_3_offset_3 >
class Periodic_3_triangulation_remove_traits_3
  : public P3DTTraits_::K
{
public:
  typedef P3DTTraits_                                           PT;
  typedef typename P3DTTraits_::K                               Base;
  typedef Off_                                                  Offset;
  typedef Periodic_3_triangulation_remove_traits_3<PT, Offset>  Self;

  typedef typename PT::RT                RT;
  typedef typename PT::FT                FT;
  typedef typename PT::Point_3           Bare_point;
  typedef std::pair<Bare_point,Offset>       Point_3;
  typedef typename PT::Iso_cuboid_3      Iso_cuboid_3;

  Periodic_3_triangulation_remove_traits_3(const Iso_cuboid_3& domain)
    : _pt()
  {
    _pt.set_domain(domain);
  }

  // Construct point
  typedef Construct_point_3_remove_traits<Self> Construct_point_3;

  // Triangulation traits
  typedef Point_offset_adaptor<Self, typename PT::Compare_xyz_3>
      Compare_xyz_3;
  typedef Point_offset_adaptor<Self, typename PT::Coplanar_orientation_3>
      Coplanar_orientation_3;
  typedef Point_offset_adaptor<Self, typename PT::Orientation_3>
      Orientation_3;

  // Delaunay Triangulation traits
  typedef Point_offset_adaptor<Self,
          typename PT::Coplanar_side_of_bounded_circle_3>
      Coplanar_side_of_bounded_circle_3;
  typedef Point_offset_adaptor<Self, typename PT::Side_of_oriented_sphere_3>
      Side_of_oriented_sphere_3;
  typedef Point_offset_adaptor<Self, typename PT::Compare_distance_3>
      Compare_distance_3;

  // Operations
  Construct_point_3
  construct_point_3_object() const {
    return Construct_point_3();
  }

  Compare_xyz_3
  compare_xyz_3_object() const {
    return Compare_xyz_3(_pt.compare_xyz_3_object());
  }
  Coplanar_orientation_3
  coplanar_orientation_3_object() const {
    return Coplanar_orientation_3(_pt.coplanar_orientation_3_object());
  }
  Orientation_3
  orientation_3_object() const {
    return Orientation_3(_pt.orientation_3_object());
  }
  Coplanar_side_of_bounded_circle_3
  coplanar_side_of_bounded_circle_3_object() const {
    return Coplanar_side_of_bounded_circle_3(_pt.coplanar_side_of_bounded_circle_3_object());
  }
  Side_of_oriented_sphere_3
  side_of_oriented_sphere_3_object() const {
    return Side_of_oriented_sphere_3(_pt.side_of_oriented_sphere_3_object());
  }
  Compare_distance_3
  compare_distance_3_object() const {
    return Compare_distance_3(_pt.compare_distance_3_object());
  }
public:
  PT _pt;
};

} //namespace CGAL

#endif // CGAL_PERIODIC_3_TRIANGULATION_REMOVE_TRAITS_3_H
