// Copyright (c) 1997, 2012, 2019 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the Licenxse, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Mariette Yvinec, Claudia Werner

#ifndef CGAL_PROJECTION_SPHERE_TRAITS_3_H
#define CGAL_PROJECTION_SPHERE_TRAITS_3_H

#include <CGAL/Delaunay_triangulation_sphere_traits_2.h>
#include <CGAL/triangulation_assertions.h>

#include <CGAL/number_utils_classes.h>
#include <CGAL/Kernel_traits.h>

namespace CGAL {

template < typename K>
class Point_with_scale
  : public K::Point_3
{
public:
  typedef typename K::FT                              FT;
  typedef typename K::Point_3                         Base_point;

  Point_with_scale() : Base_point() { }
  Point_with_scale(const Base_point& p,
                   const typename K::Point_3& sphere_center)
    : Base_point(p)
  {
    compute_scale(p-(sphere_center-ORIGIN));
  }

  FT scale() const { return _scale; }

private:
  void compute_scale(const FT x, const FT y, const FT z)
  {
    const FT tmp = x*x + y*y + z*z;
    if(tmp == 0)
      _scale = 0;
    else
      _scale = 1 / CGAL::sqrt(tmp);
  }

  void compute_scale(const Base_point& p)
  {
    return compute_scale(p.x(), p.y(), p.z());
  }

  FT _scale;
};

//adaptor for calling the Predicate_ with the points projected on the sphere
template <typename Traits, typename Predicate_>
class Functor_projection_adaptor
  : public Predicate_
{
  typedef Predicate_                                 Predicate;
  typedef Predicate_                                 Base;

  typedef typename Traits::FT                        FT;
  typedef typename Traits::Point_2                   Point; // that's Projected_point<>
  typedef typename Traits::Base_point                Base_point;

public:
  typedef typename Predicate::result_type            result_type;

  Functor_projection_adaptor(const Predicate& p, const FT r) : Base(p), _radius(r) { }

  result_type operator()(const Point& p0, const Point& p1)
  { return Base::operator()(project(p0), project(p1)); }

  result_type operator()(const Point& p0, const Point& p1, const Point& p2)
  { return Base::operator()(project(p0), project(p1), project(p2)); }

  result_type operator ()(const Point& p0, const Point& p1, const Point& p2, const Point& p3)
  { return Base::operator()(project(p0), project(p1), project(p2), project(p3)); }

  result_type operator()(const Point& p0, const Point& p1, const Point& p2, const Point& p3, const Point& p4)
  { return Base::operator()(project(p0), project(p1), project(p2), project(p3), project(p4)); }

private:
  Base_point project(const Point& p)
  {
    const FT scale = _radius * p.scale();
    return Base_point(scale*p.x(), scale*p.y(), scale*p.z());
  }

  const FT _radius;
};

template <typename K>
class Projection_sphere_traits_3
  : public Delaunay_triangulation_sphere_traits_2<K>
{
protected:
  typedef Delaunay_triangulation_sphere_traits_2<K>                          Base;
  typedef Projection_sphere_traits_3<K>                                      Self;

public:
  typedef typename K::FT                                                     FT;
  typedef Point_with_scale<K>                                                Point_2;
  typedef typename K::Point_3                                                Base_point;

  typedef Functor_projection_adaptor<Self, typename Base::Power_test_2>      Power_test_2;
  typedef Functor_projection_adaptor<Self, typename Base::Orientation_2>     Orientation_2;
  typedef Functor_projection_adaptor<Self, typename Base::Coradial_sphere_2> Coradial_sphere_2;
  typedef Functor_projection_adaptor<Self, typename Base::Inside_cone_2>     Inside_cone_2;

  typedef Functor_projection_adaptor<Self, typename Base::Orientation_1>     Orientation_1;
  typedef Functor_projection_adaptor<Self, typename Base::Compute_squared_distance_2> Compute_squared_distance_2;
  typedef Functor_projection_adaptor<Self, typename Base::Compare_xyz_3>     Compare_xyz_3;

  Projection_sphere_traits_3(const Base_point& sphere = CGAL::ORIGIN,
                             const FT radius = 1,
                             const K& k = K())
    : Base(sphere, radius, k)
  { }

  void set_radius(double radius) { _radius = radius; }

  Orientation_2
  orientation_2_object() const
  { return Orientation_2(Base::orientation_2_object(), _radius); }

  Orientation_1
  orientation_1_object() const
  { return Orientation_1(Base::orientation_1_object(), _radius); }

  Power_test_2
  power_test_2_object() const
  { return Power_test_2(Base::power_test_2_object(), _radius); }

  Coradial_sphere_2
  coradial_sphere_2_object() const
  { return Coradial_sphere_2(Base::coradial_sphere_2_object(), _radius); }

  Inside_cone_2
  inside_cone_2_object() const
  { return Inside_cone_2(Base::inside_cone_2_object(), _radius); }

  Compute_squared_distance_2
  compute_squared_distance_3_object() const
  { return Compute_squared_distance_2(Base::compute_squared_distance_2_object(), _radius); }

  Compare_xyz_3
  compare_xyz_3_object() const
  { return Compare_xyz_3(Base::compare_xyz_3_object(), _radius); }

  // @fixme rename that?
  struct Construct_projected_point_3
    : public std::unary_function<Base_point, Point_2>
  {
    Construct_projected_point_3(const Base_point& sc) : sphere_center(sc) { }

    Point_2 operator()(const Base_point& pt) const { return Point_2(pt, sphere_center); }

  private:
    const Base_point& sphere_center;
  };

  Construct_projected_point_3
  construct_projected_point_3_object() const
  { return Construct_projected_point_3(_sphere); }

protected :
  FT _radius;
  Base_point _sphere;
};

} // namespace CGAL

#endif // CGAL_PROJECTION_SPHERE_TRAITS_3_H
