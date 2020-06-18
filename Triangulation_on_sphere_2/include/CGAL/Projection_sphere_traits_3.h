// Copyright (c) 1997, 2012, 2019 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mariette Yvinec, Claudia Werner

#ifndef CGAL_TRIANGULATION_ON_SPHERE_PROJECTION_SPHERE_TRAITS_3_H
#define CGAL_TRIANGULATION_ON_SPHERE_PROJECTION_SPHERE_TRAITS_3_H

#include <CGAL/license/Triangulation_on_sphere_2.h>

#include <CGAL/Delaunay_triangulation_sphere_traits_2.h>

#include <CGAL/triangulation_assertions.h>
#include <CGAL/number_utils_classes.h>

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
      _scale = FT(0);
    else
      _scale = FT(1) / CGAL::approximate_sqrt(tmp);
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
  typedef typename Traits::Point_on_sphere_2         Point;
  typedef typename Traits::Point_3                   Point_3;

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
  typedef Delaunay_triangulation_sphere_traits_2<K>                                   Base;
  typedef Projection_sphere_traits_3<K>                                               Self;

public:
  typedef typename K::FT                                                              FT;
  typedef Point_with_scale<K>                                                         Point_on_sphere_2;
  typedef typename K::Point_3                                                         Point_3;

  typedef Functor_projection_adaptor<Self, typename Base::Construct_circumcenter_on_sphere_2>  Construct_circumcenter_on_sphere_2;
  typedef Functor_projection_adaptor<Self, typename Base::Equal_on_sphere_2>                   Equal_on_sphere_2;
  typedef Functor_projection_adaptor<Self, typename Base::Inside_cone_2>                       Inside_cone_2;
  typedef Functor_projection_adaptor<Self, typename Base::Orientation_on_sphere_2>             Orientation_on_sphere_2;
  typedef Functor_projection_adaptor<Self, typename Base::Power_test_2>                        Power_test_2;

  Projection_sphere_traits_3(const Point_3& sphere = CGAL::ORIGIN,
                             const FT radius = 1,
                             const K& k = K())
    : Base(sphere, radius, k)
  { }

public:
  Construct_circumcenter_on_sphere_2
  construct_circumcenter_on_sphere_2_object() const
  { return Construct_circumcenter_on_sphere_2(Base::construct_circumcenter_on_sphere_2(), Base::_radius); }

  Equal_on_sphere_2
  equal_on_sphere_2_object() const
  { return Equal_on_sphere_2(Base::equal_on_sphere_2(), Base::_radius); }

  Inside_cone_2
  inside_cone_2_object() const
  { return Inside_cone_2(Base::inside_cone_2_object(), Base::_radius); }

  Orientation_on_sphere_2
  orientation_on_sphere_2() const
  { return Orientation_2(Base::orientation_on_sphere_2_object(), Base::_radius); }

  Power_test_2
  power_test_2_object() const
  { return Power_test_2(Base::power_test_2_object(), Base::_radius); }

public:
  // @fixme rename that?
  struct Construct_projected_point_3
    : public std::unary_function<Point_3, Point_on_sphere_2>
  {
    Construct_projected_point_3(const Point_3& sc) : sphere_center(sc) { }

    Point_on_sphere_2 operator()(const Point_3& pt) const { return Point_on_sphere_2(pt, sphere_center); }

  private:
    const Point_3& sphere_center;
  };

  Construct_projected_point_3
  construct_projected_point_3_object() const
  { return Construct_projected_point_3(_center); }
};

} // namespace CGAL

#endif // CGAL_TRIANGULATION_ON_SPHERE_PROJECTION_SPHERE_TRAITS_3_H
