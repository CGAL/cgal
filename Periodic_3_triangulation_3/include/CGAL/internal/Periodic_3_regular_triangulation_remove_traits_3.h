// Copyright (c) 1999-2004,2006-2009,2013-2015,2017  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Monique Teillaud <Monique.Teillaud@inria.fr>
//                 Aymeric Pelle <Aymeric.Pelle@sophia.inria.fr>

#ifndef CGAL_PERIODIC_3_REGULAR_TRIANGULATION_REMOVE_TRAITS_3_H
#define CGAL_PERIODIC_3_REGULAR_TRIANGULATION_REMOVE_TRAITS_3_H

#include <CGAL/license/Periodic_3_triangulation_3.h>

#include <CGAL/basic.h>
#include <CGAL/Periodic_3_offset_3.h>

#include <CGAL/internal/Periodic_3_Delaunay_triangulation_remove_traits_3.h>

#include <utility>

namespace CGAL {

// Triangulation_3 has calls to Construct_point_3 to handle weighted and bare points.
// The default inherited Construct_point_3 inherited by Periodic_3_triangulation_remove_traits_3
// must be overwritten by a construction Construct_point_3 that offers:
// - pair<K::Point_3, offset> --> pair<K::Point_3, offset> (identity)
// - pair<K::Weighted_point_3, offset> --> pair<K::Point_3, offset>
template<class Gt_,
         typename Construct_point_3_base_>
class Construct_point_from_weighted_pair_3
  : public Construct_point_from_pair_3<Gt_, Construct_point_3_base_>
{
  typedef Construct_point_from_pair_3<Gt_, Construct_point_3_base_>  Base;

  typedef Gt_                                                        Geom_traits;

  // `Traits::Point_3` is actually a `std::pair<Point_3, Offset>`
  // `Traits::Weighted_point_3` is actually a `std::pair<Weighted_point_3, Offset>`
  typedef typename Geom_traits::Point_3                              Point_3;
  typedef typename Geom_traits::Weighted_point_3                     Weighted_point_3;

public:
  Construct_point_from_weighted_pair_3(const Construct_point_3_base_& cp) : Base(cp) { }

  using Base::operator(); // for K::Weighted_point_3 to Point_3

  template<typename F>
  struct result : Base::template result<F> {};

  template<typename F>
  struct result<F(Weighted_point_3)> {
    typedef Point_3 type;
  };

  Point_3 operator()(const Weighted_point_3& wp) const {
    return std::make_pair(operator()(wp.first), wp.second /* offset */);
  }
};

template<class Traits_, class Functor_>
class Functor_with_weighted_point_offset_pair_adaptor
  : public Functor_
{
  typedef Traits_                                 Traits;
  typedef Functor_                                Functor;

  // `Traits::Point_3` is actually a `std::pair<Point_3, Offset>`
  // `Traits::Weighted_point_3` is actually a `std::pair<Weighted_point_3, Offset>`
  typedef typename Traits::Point_3                Point_3;
  typedef typename Traits::Weighted_point_3       Weighted_point_3;

public:
  typedef typename Functor::result_type result_type;

  Functor_with_weighted_point_offset_pair_adaptor (const Functor & functor)
    : Functor_(functor)
  { }

public:
  using Functor::operator();

  result_type operator() (const Point_3& p0, const Weighted_point_3& p1,
                          const Weighted_point_3& p2) const
  {
    return operator()(p0.first, p1.first, p2.first,
                      p0.second, p1.second, p2.second);
  }

  // bare points
  result_type operator()(const Point_3& p0, const Point_3& p1) const {
    return operator()(p0.first, p1.first,
                      p0.second, p1.second);
  }
  result_type operator()(const Point_3& p0, const Point_3& p1,
                         const Point_3& p2) const {
    return operator()(p0.first, p1.first, p2.first,
                      p0.second, p1.second, p2.second);
  }
  result_type operator()(const Point_3& p0, const Point_3& p1,
                         const Point_3& p2, const Point_3& p3) const {
    return operator()(p0.first, p1.first, p2.first, p3.first,
                      p0.second, p1.second, p2.second, p3.second);
  }
  result_type operator()(const Point_3& p0, const Point_3& p1,
                         const Point_3& p2, const Point_3& p3, const Point_3& p4) const {
    return operator()(p0.first, p1.first, p2.first, p3.first, p4.first,
                      p0.second, p1.second, p2.second, p3.second, p4.second);
  }

  // weighted points
  result_type operator() (const Weighted_point_3& p0, const Weighted_point_3& p1) const
  {
    return operator()(p0.first, p1.first, p0.second, p1.second);
  }
  result_type operator() (const Weighted_point_3& p0, const Weighted_point_3& p1,
                          const Weighted_point_3& p2) const
  {
    return operator()(p0.first, p1.first, p2.first,
                      p0.second, p1.second, p2.second);
  }
  result_type operator() (const Weighted_point_3& p0, const Weighted_point_3& p1,
                          const Weighted_point_3& p2, const Weighted_point_3& p3) const
  {
    return operator()(p0.first, p1.first, p2.first, p3.first,
                      p0.second, p1.second, p2.second, p3.second);
  }
  result_type operator() (const Weighted_point_3& p0, const Weighted_point_3& p1,
                          const Weighted_point_3& p2, const Weighted_point_3& p3,
                          const Weighted_point_3& p4) const
  {
    return operator()(p0.first, p1.first, p2.first, p3.first, p4.first,
                      p0.second, p1.second, p2.second, p3.second, p4.second);
  }
};

template<class Gt_, class Off_ = typename CGAL::Periodic_3_offset_3>
class Periodic_3_regular_triangulation_remove_traits_3
  : public Gt_
{
  typedef Periodic_3_regular_triangulation_remove_traits_3<Gt_, Off_>    Self;
  typedef Gt_                                                            Base;

public:
  typedef Gt_                                                            Geom_traits;
  typedef Off_                                                           Offset;

  typedef typename Geom_traits::RT                                       RT;
  typedef typename Geom_traits::FT                                       FT;
  typedef std::pair<typename Geom_traits::Point_3, Offset>               Point_3;
  typedef std::pair<typename Geom_traits::Weighted_point_3, Offset>      Weighted_point_3;

  // not allowing a default value for `gt` because we need to have
  // an initialized domain in `gt`
  Periodic_3_regular_triangulation_remove_traits_3(const Geom_traits& gt) : Base(gt) { }

  // Construct point
  typedef Construct_point_from_weighted_pair_3<Self, typename Geom_traits::Construct_point_3>
      Construct_point_3;

  // Triangulation traits
  typedef Functor_with_weighted_point_offset_pair_adaptor<Self, typename Geom_traits::Compare_xyz_3>
      Compare_xyz_3;
  typedef Functor_with_weighted_point_offset_pair_adaptor<Self, typename Geom_traits::Coplanar_orientation_3>
      Coplanar_orientation_3;
  typedef Functor_with_weighted_point_offset_pair_adaptor<Self, typename Geom_traits::Orientation_3>
      Orientation_3;

  // Regular Triangulation traits
  typedef Functor_with_weighted_point_offset_pair_adaptor<Self, typename Geom_traits::Power_side_of_oriented_power_sphere_3>
      Power_side_of_oriented_power_sphere_3;
  typedef Functor_with_weighted_point_offset_pair_adaptor<Self, typename Geom_traits::Compare_power_distance_3>
      Compare_power_distance_3;
  typedef Functor_with_weighted_point_offset_pair_adaptor<Self, typename Geom_traits::Power_side_of_bounded_power_sphere_3>
      Power_side_of_bounded_power_sphere_3;
  typedef Functor_with_weighted_point_offset_pair_adaptor<Self, typename Geom_traits::Construct_weighted_circumcenter_3>
      Construct_weighted_circumcenter_3;
  typedef Functor_with_weighted_point_offset_pair_adaptor<Self, typename Geom_traits::Construct_circumcenter_3>
      Construct_circumcenter_3;
  typedef Functor_with_weighted_point_offset_pair_adaptor<Self, typename Geom_traits::Compute_squared_radius_smallest_orthogonal_sphere_3>
      Compute_squared_radius_smallest_orthogonal_sphere_3;
  typedef Functor_with_weighted_point_offset_pair_adaptor<Self, typename Geom_traits::Compute_power_product_3>
      Compute_power_product_3;
  typedef Functor_with_weighted_point_offset_pair_adaptor<Self, typename Geom_traits::Compute_power_distance_to_power_sphere_3>
      Compute_power_distance_to_power_sphere_3;
  typedef Functor_with_weighted_point_offset_pair_adaptor<Self, typename Geom_traits::Compare_weighted_squared_radius_3>
      Compare_weighted_squared_radius_3;

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
  Power_side_of_oriented_power_sphere_3
  power_side_of_oriented_power_sphere_3_object() const {
    return Power_side_of_oriented_power_sphere_3(this->Base::power_side_of_oriented_power_sphere_3_object());
  }
  Compare_power_distance_3
  compare_power_distance_3_object() const {
    return Compare_power_distance_3(this->Base::compare_power_distance_3_object());
  }
  Power_side_of_bounded_power_sphere_3
  power_side_of_bounded_power_sphere_3_object() const {
    return Power_side_of_bounded_power_sphere_3(this->Base::power_side_of_bounded_power_sphere_3_object());
  }
  Construct_weighted_circumcenter_3
  construct_weighted_circumcenter_3_object() const {
    return Construct_weighted_circumcenter_3(this->Base::construct_weighted_circumcenter_3_object());
  }
  Construct_circumcenter_3
  construct_circumcenter_3_object() const {
    return Construct_circumcenter_3(this->Base::construct_circumcenter_3_object());
  }
  Compute_squared_radius_smallest_orthogonal_sphere_3
  compute_squared_radius_smallest_orthogonal_sphere_3_object() const {
    return Compute_squared_radius_smallest_orthogonal_sphere_3(this->Base::compute_squared_radius_smallest_orthogonal_sphere_3_object());
  }
  Compute_power_product_3
  compute_power_product_3_object() const {
    return Compute_power_product_3(this->Base::compute_power_product_3_object());
  }
  Compute_power_distance_to_power_sphere_3
  compute_power_distance_to_power_sphere_3_object() const {
    return Compute_power_distance_to_power_sphere_3(this->Base::compute_power_distance_to_power_sphere_3_object());
  }
  Compare_weighted_squared_radius_3
  compare_weighted_squared_radius_3_object() const {
    return Compare_weighted_squared_radius_3(this->Base::compare_weighted_squared_radius_3_object());
  }
};

} //namespace CGAL

#endif // CGAL_PERIODIC_3_REGULAR_TRIANGULATION_REMOVE_TRAITS_3_H
