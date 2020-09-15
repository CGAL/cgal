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

#ifndef CGAL_PERIODIC_3_DELAUNAY_TRIANGULATION_REMOVE_TRAITS_3_H
#define CGAL_PERIODIC_3_DELAUNAY_TRIANGULATION_REMOVE_TRAITS_3_H

#include <CGAL/license/Periodic_3_triangulation_3.h>

#include <CGAL/basic.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Periodic_3_offset_3.h>

namespace CGAL {

// Triangulation_3 uses Construct_point_3 to handle weighted and bare points.
// The default Construct_point_3 inherited by Periodic_3_triangulation_remove_traits_3
// must be overwritten by a custom Construct_point_3 that offers:
// - pair<K::Point_3, offset> --> pair<K::Point_3, offset> (identity)
template<class Gt_, typename Construct_point_3_base_>
class Construct_point_from_pair_3
  : public Construct_point_3_base_
{
  typedef Construct_point_3_base_        Base;
  typedef Gt_                            Geom_traits;

  // Gt::Point_3 is actually a pair <K::Point_3, offset>
  typedef typename Geom_traits::Point_3  Point_3;

public:
  Construct_point_from_pair_3(const Base& cp) : Base(cp) { }

  using Base::operator();

  template<typename F>
  struct result : Base::template result<F> {};

  template<typename F>
  struct result<F(Point_3)> {
    typedef const Point_3& type;
  };

  const Point_3& operator()(const Point_3& p) const { return p; }
};

template < class Traits_, class Functor_ >
class Functor_with_point_offset_pair_adaptor
  : public Functor_
{
  typedef Traits_                        Traits;
  typedef Functor_                       Functor;

  // `Traits::Point_3` is actually a `std::pair<Point_3, Offset>`
  typedef typename Traits::Point_3       Point;

public:
  typedef typename Functor::result_type result_type;

  Functor_with_point_offset_pair_adaptor(const Functor & functor) : Functor_(functor) { }

public:
  using Functor::operator();

  result_type operator()(const Point& p0, const Point& p1) const {
    return operator()(p0.first, p1.first,
                      p0.second, p1.second);
  }
  result_type operator()(const Point& p0, const Point& p1,
                         const Point& p2) const {
    return operator()(p0.first, p1.first, p2.first,
                      p0.second, p1.second, p2.second);
  }
  result_type operator()(const Point& p0, const Point& p1,
                         const Point& p2, const Point& p3) const {
    return operator()(p0.first, p1.first, p2.first, p3.first,
                      p0.second, p1.second, p2.second, p3.second);
  }
  result_type operator()(const Point& p0, const Point& p1,
                         const Point& p2, const Point& p3, const Point& p4) const {
    return operator()(p0.first, p1.first, p2.first, p3.first, p4.first,
                      p0.second, p1.second, p2.second, p3.second, p4.second);
  }
};

template < class Gt_, class Off_ = typename CGAL::Periodic_3_offset_3 >
class Periodic_3_Delaunay_triangulation_remove_traits_3
  : public Gt_
{
  typedef Periodic_3_Delaunay_triangulation_remove_traits_3<Gt_, Off_>  Self;
  typedef Gt_                                                           Base;

public:
  typedef Gt_                                                    Geom_traits;
  typedef Off_                                                   Offset;

  typedef typename Geom_traits::RT                               RT;
  typedef typename Geom_traits::FT                               FT;
  typedef std::pair<typename Geom_traits::Point_3, Offset>       Point_3;

  // not allowing a default value for `gt` because we need to have
  // an initialized domain in `gt`
  Periodic_3_Delaunay_triangulation_remove_traits_3(const Geom_traits& gt) : Base(gt) { }

  // Construct point
  typedef Construct_point_from_pair_3<Self, typename Geom_traits::Construct_point_3> Construct_point_3;

  // Triangulation predicates
  typedef Functor_with_point_offset_pair_adaptor<Self, typename Geom_traits::Compare_xyz_3>
      Compare_xyz_3;
  typedef Functor_with_point_offset_pair_adaptor<Self, typename Geom_traits::Orientation_3>
      Orientation_3;

  // Delaunay Triangulation predicates
  typedef Functor_with_point_offset_pair_adaptor<Self, typename Geom_traits::Compare_distance_3>
      Compare_distance_3;
  typedef Functor_with_point_offset_pair_adaptor<Self, typename Geom_traits::Side_of_oriented_sphere_3>
      Side_of_oriented_sphere_3;

  // Degenerate dimension predicates
  typedef Functor_with_point_offset_pair_adaptor<Self, typename Geom_traits::Coplanar_orientation_3>
      Coplanar_orientation_3;
  typedef Functor_with_point_offset_pair_adaptor<Self, typename Geom_traits::Coplanar_side_of_bounded_circle_3>
      Coplanar_side_of_bounded_circle_3;

  // Operations
  Construct_point_3
  construct_point_3_object() const {
    return Construct_point_3(this->Base::construct_point_3_object());
  }

  Compare_xyz_3
  compare_xyz_3_object() const {
    return Compare_xyz_3(this->Base::compare_xyz_3_object());
  }
  Coplanar_orientation_3
  coplanar_orientation_3_object() const {
    return Coplanar_orientation_3(this->Base::coplanar_orientation_3_object());
  }
  Orientation_3
  orientation_3_object() const {
    return Orientation_3(this->Base::orientation_3_object());
  }
  Coplanar_side_of_bounded_circle_3
  coplanar_side_of_bounded_circle_3_object() const {
    return Coplanar_side_of_bounded_circle_3(this->Base::coplanar_side_of_bounded_circle_3_object());
  }
  Side_of_oriented_sphere_3
  side_of_oriented_sphere_3_object() const {
    return Side_of_oriented_sphere_3(this->Base::side_of_oriented_sphere_3_object());
  }
  Compare_distance_3
  compare_distance_3_object() const {
    return Compare_distance_3(this->Base::compare_distance_3_object());
  }
};

} //namespace CGAL

#endif // CGAL_PERIODIC_3_DELAUNAY_TRIANGULATION_REMOVE_TRAITS_3_H
