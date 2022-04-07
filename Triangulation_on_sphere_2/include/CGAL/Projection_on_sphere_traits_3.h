// Copyright (c) 1997, 2012, 2019 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mariette Yvinec,
//                 Claudia Werner
//                 Mael Rouxel-Labbé

#ifndef CGAL_TRIANGULATION_ON_SPHERE_PROJECTION_SPHERE_TRAITS_3_H
#define CGAL_TRIANGULATION_ON_SPHERE_PROJECTION_SPHERE_TRAITS_3_H

#include <CGAL/license/Triangulation_on_sphere_2.h>

#include <CGAL/Delaunay_triangulation_on_sphere_traits_2.h>

#include <CGAL/triangulation_assertions.h>
#include <CGAL/number_utils_classes.h>

namespace CGAL {
namespace internal {

template <typename LK>
class Point_with_scale
  : public LK::Point_3
{
public:
  typedef typename LK::FT                                  FT;
  typedef typename LK::Point_3                             Base_point;
  typedef typename LK::Vector_3                            Vector_3;

private:
  mutable bool cached;
  mutable Base_point proj_pt;

public:
  Point_with_scale() : Base_point() { } // vertex base wants a default constructor
  Point_with_scale(const Base_point& p) : Base_point(p), cached(false) { }

  const Base_point& get_projection(const Base_point& center,
                                   const FT radius,
                                   const LK& lk) const
  {
    if(!cached)
      compute_projection(center, radius, lk);

    return proj_pt;
  }

  const Base_point& get_projection(const Base_point& center,
                                   const FT radius) const
  {
    return get_projection(center, radius, LK());
  }

private:
  void compute_projection(const Base_point& center,
                          const FT radius,
                          const LK& lk) const
  {
    CGAL_precondition(static_cast<const Base_point&>(*this) != center);

    Vector_3 v = lk.construct_vector_3_object()(center, static_cast<const Base_point&>(*this));
    FT n = lk.compute_squared_length_3_object()(v);
    CGAL_assertion(n > 0);
    n = CGAL::sqrt(n);
    v = lk.construct_scaled_vector_3_object()(v, radius / n);

    proj_pt = lk.construct_translated_point_3_object()(center, v);
    cached = true;
  }
};

template <typename P3, typename PoS2>
struct Construct_point_with_scale
  : public CGAL::cpp98::unary_function<P3, PoS2>
{
  PoS2 operator()(const P3& pt) const { return PoS2(pt); }
};

// Adaptor for calling the functors with the points projected on the sphere
template <typename LK, typename Functor_>
class Functor_projection_adaptor
  : public Functor_
{
  typedef Functor_                                         Functor;
  typedef Functor_                                         Base;

  typedef typename LK::FT                                  FT;
  typedef typename LK::Point_3                             Point_3;
  typedef internal::Point_with_scale<LK>                   Point;

public:
  Functor_projection_adaptor(const Functor& f, const Point_3& center, const FT radius)
    : Base(f), _c(center), _r(radius)
  { }

public:
  using Base::operator();

  decltype(auto) operator()(const Point& p0, const Point& p1)
  {
    return Base::operator()(p0.get_projection(_c, _r), p1.get_projection(_c, _r));
  }

  decltype(auto) operator()(const Point& p0, const Point& p1, const Point& p2)
  {
    return Base::operator()(p0.get_projection(_c, _r), p1.get_projection(_c, _r),
                            p2.get_projection(_c, _r));
  }

  decltype(auto) operator ()(const Point& p0, const Point& p1, const Point& p2, const Point& p3)
  {
    return Base::operator()(p0.get_projection(_c, _r), p1.get_projection(_c, _r),
                            p2.get_projection(_c, _r), p3.get_projection(_c, _r));
  }

  decltype(auto) operator()(const Point& p0, const Point& p1, const Point& p2,
                            const Point& p3, const Point& p4)
  {
    return Base::operator()(p0.get_projection(_c, _r), p1.get_projection(_c, _r),
                            p2.get_projection(_c, _r), p3.get_projection(_c, _r),
                            p4.get_projection(_c, _r));
  }

private:
  const Point_3& _c;
  const FT _r;
};

} // namespace internal

template <typename LK>
class Projection_on_sphere_traits_3
  : public Delaunay_triangulation_on_sphere_traits_2<LK>
{
protected:
  typedef Delaunay_triangulation_on_sphere_traits_2<LK>                        Base;
  typedef Projection_on_sphere_traits_3<LK>                                    Self;

public:
  typedef typename LK::FT                                                      FT;
  typedef typename LK::Point_3                                                 Point_3;
  typedef internal::Point_with_scale<LK>                                       Point_on_sphere_2;

  // predicates
  typedef internal::Functor_projection_adaptor<LK, typename Base::Collinear_are_strictly_ordered_on_great_circle_2> Collinear_are_strictly_ordered_on_great_circle_2;
  typedef internal::Functor_projection_adaptor<LK, typename Base::Compare_on_sphere_2> Compare_on_sphere_2;
  typedef internal::Functor_projection_adaptor<LK, typename Base::Equal_on_sphere_2> Equal_on_sphere_2;
  typedef internal::Functor_projection_adaptor<LK, typename Base::Orientation_on_sphere_2> Orientation_on_sphere_2;
  typedef internal::Functor_projection_adaptor<LK, typename Base::Side_of_oriented_circle_on_sphere_2> Side_of_oriented_circle_on_sphere_2;

  // constructions
  typedef internal::Construct_point_with_scale<Point_3, Point_on_sphere_2> Construct_point_on_sphere_2;
  typedef internal::Functor_projection_adaptor<LK, typename Base::Construct_point_3> Construct_point_3;
  typedef internal::Functor_projection_adaptor<LK, typename Base::Construct_segment_3> Construct_segment_3;
  typedef internal::Functor_projection_adaptor<LK, typename Base::Construct_triangle_3> Construct_triangle_3;

  typedef internal::Functor_projection_adaptor<LK, typename Base::Construct_arc_on_sphere_2> Construct_arc_on_sphere_2;
  typedef internal::Functor_projection_adaptor<LK, typename Base::Construct_circumcenter_on_sphere_2> Construct_circumcenter_on_sphere_2;

public:
  Projection_on_sphere_traits_3(const Point_3& sphere = CGAL::ORIGIN,
                                const FT radius = 1,
                                const LK& lk = LK())
    : Base(sphere, radius, lk)
  { }

  // predicates
public:
  Collinear_are_strictly_ordered_on_great_circle_2
  collinear_are_strictly_ordered_on_great_circle_2_object() const
  { return Collinear_are_strictly_ordered_on_great_circle_2(
             Base::collinear_are_strictly_ordered_on_great_circle_2_object(),
             Base::center(), Base::radius()); }

  Compare_on_sphere_2
  compare_on_sphere_2_object() const
  { return Compare_on_sphere_2(Base::compare_on_sphere_2_object(),
                               Base::center(), Base::radius()); }

  Equal_on_sphere_2
  equal_on_sphere_2_object() const
  { return Equal_on_sphere_2(Base::equal_on_sphere_2_object(),
                             Base::center(), Base::radius()); }

  Orientation_on_sphere_2
  orientation_on_sphere_2_object() const
  { return Orientation_on_sphere_2(Base::orientation_on_sphere_2_object(),
                                   Base::center(), Base::radius()); }

  Side_of_oriented_circle_on_sphere_2
  side_of_oriented_circle_on_sphere_2_object() const
  { return Side_of_oriented_circle_on_sphere_2(Base::side_of_oriented_circle_on_sphere_2_object(),
                                               Base::center(), Base::radius()); }

  // constructions
public:
  Construct_point_on_sphere_2
  construct_point_on_sphere_2_object() const
  { return Construct_point_on_sphere_2(); }

  Construct_arc_on_sphere_2
  construct_arc_on_sphere_2_object() const
  { return Construct_arc_on_sphere_2(Base::construct_arc_on_sphere_2_object(),
                                     Base::center(), Base::radius()); }

  Construct_circumcenter_on_sphere_2
  construct_circumcenter_on_sphere_2_object() const
  { return Construct_circumcenter_on_sphere_2(Base::construct_circumcenter_on_sphere_2_object(),
                                              Base::center(), Base::radius()); }

  Construct_point_3
  construct_point_3_object() const
  { return Construct_point_3(Base::construct_point_3_object(), Base::center(), Base::radius()); }

public:
  bool is_on_sphere(const Point_on_sphere_2& /*p*/) const
  {
    return true;
  }

  bool are_points_too_close(const Point_on_sphere_2& p, const Point_on_sphere_2& q) const
  {
    return Base::are_points_too_close(p.get_projection(Base::center(), Base::radius(), Base::lk()),
                                      q.get_projection(Base::center(), Base::radius(), Base::lk()));
  }
};

} // namespace CGAL

#endif // CGAL_TRIANGULATION_ON_SPHERE_PROJECTION_SPHERE_TRAITS_3_H
