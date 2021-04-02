// Copyright (c) 1997, 2012, 2019 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_TRIANGULATION_ON_SPHERE_GEOGRAPHICAL_COORDINATES_TRAITS_2_H
#define CGAL_TRIANGULATION_ON_SPHERE_GEOGRAPHICAL_COORDINATES_TRAITS_2_H

#include <CGAL/license/Triangulation_on_sphere_2.h>

#include <CGAL/Delaunay_triangulation_on_sphere_traits_2.h>

#include <CGAL/triangulation_assertions.h>
#include <CGAL/number_utils_classes.h>

#include <cmath>
#include <fstream>
#include <iostream>

// @todo
// Note: this is currently undocumented because it uses a conversion to
//       the Cartesian domain R^3 to perform predicates and constructions,
//       losing the benefit of the exact representation on the way.
//       It can be useful on its own (especially since it contains code
//       to handle generic convex domains), but hidden points are possible
//       and are not handled properly.

namespace CGAL {

template <typename LK>
struct Geographical_coordinates
{
  typedef typename LK::FT                           FT;
  typedef FT                                        Latitude;
  typedef FT                                        Longitude;

  Geographical_coordinates() : _la(360), _lo(360) { }
  Geographical_coordinates(const FT la, const FT lo) : _la(la), _lo(lo) { }

  Latitude latitude() const { return _la; }
  Longitude longitude() const { return _lo; }

  friend std::ostream& operator<<(std::ostream& os, const Geographical_coordinates& l)
  {
    os << l._la << " " << l._lo;
    return os;
  }

private:
  FT _la;
  FT _lo;
};

namespace internal {

template <typename LK>
struct Construct_geographical_coordinates
  : public CGAL::cpp98::unary_function<typename LK::Point_3, CGAL::Geographical_coordinates<LK> >
{
  typedef typename LK::FT                                            FT;
  typedef typename LK::Point_3                                       Point_3;
  typedef CGAL::Geographical_coordinates<LK>                         result_type;

  Construct_geographical_coordinates(const Point_3& center, const FT radius)
    : _center(center), _radius(radius)
  { }

  result_type operator()(const Point_3& pt) const
  {
    CGAL_assertion(pt != _center);

    const FT la = 180. * std::asin((pt.z() - _center.z()) / _radius) / CGAL_PI;
    const FT lo = 180. * std::atan2((pt.y() - _center.y()), (pt.x() - _center.x())) / CGAL_PI;

    return result_type(la, lo);
  }

private:
  const Point_3& _center;
  const FT _radius;
};

template <typename LK>
struct Construct_Cartesian_coordinates
  : public CGAL::cpp98::unary_function<typename LK::Point_3, CGAL::Geographical_coordinates<LK> >
{
  typedef typename LK::FT                                            FT;
  typedef typename LK::Point_3                                       Point_3;
  typedef CGAL::Geographical_coordinates<LK>                         Coordinates;

  Construct_Cartesian_coordinates(const Point_3& center, const FT radius)
    : _center(center), _radius(radius)
  { }

  Point_3 operator()(const Coordinates& l) const
  {
#if 1
    const FT rla = CGAL_PI * l.latitude() / 180.;
    const FT rlo = CGAL_PI * l.longitude() / 180.;
    const FT x = _center.x() + _radius * std::cos(rla) * std::cos(rlo);
    const FT y = _center.y() + _radius * std::cos(rla) * std::sin(rlo);
    const FT z = _center.z() + _radius * std::sin(rla);
#else
    const FT rla = CGAL_PI * l.latitude() / 180.;
    const FT rlo = CGAL_PI * l.longitude() / 180.;

    // Quoting wikipedia:
    // - The Earth's equatorial radius a, or semi-major axis, is the distance from its center
    // to the equator and equals 6,378.1370 km (3,963.1906 mi)
    // - The Earth's polar radius b, or semi-minor axis, is the distance from its center to the North
    // and South Poles, and equals 6,356.7523 km (3,949.9028 mi).
    const FT equ_radius = 6378.1370;
    const FT pol_radius = 6356.7523;
    FT radius_ratio = equ_radius / pol_radius;

    const FT a2 = square(_radius * radius_ratio);
    const FT b2 = square(_radius);

    // Using https://en.wikipedia.org/wiki/Geographic_coordinate_conversion#From_geodetic_to_ECEF_coordinates
    const FT N = a2 / sqrt(a2 * square(std::cos(rla)) + b2 * square(std::sin(rla)));
    const FT h = 0;
    const FT x = _center.x() + (N + h) * std::cos(rla) * std::cos(rlo);
    const FT y = _center.y() + (N + h) * std::cos(rla) * std::sin(rlo);
    const FT z = _center.z() + (b2 * N / a2 + h) * std::sin(rla);
#endif

    return Point_3(x, y, z);
  }

private:
  const Point_3& _center;
  const FT _radius;
};

template <typename LK, typename Base_functor_>
struct Construct_circumcenter_on_geo_sphere_2
{
public:
  typedef typename LK::FT                                            FT;
  typedef typename LK::Point_3                                       Point_3;
  typedef typename LK::Vector_3                                      Vector_3;
  typedef Geographical_coordinates<LK>                               Point_on_sphere_2;
  typedef Geographical_coordinates<LK>                               result_type;

  Construct_circumcenter_on_geo_sphere_2( const Base_functor_ f, const Point_3& center, const FT radius)
    : _center(center), _radius(radius), _f(f)
  { }

  result_type operator()(const Point_on_sphere_2& p, const Point_on_sphere_2& q, const Point_on_sphere_2& r) const
  {
    Construct_Cartesian_coordinates<LK> to_P3(_center, _radius);
    Construct_geographical_coordinates<LK> to_PoS2(_center, _radius);

    return to_PoS2(_f(to_P3(p), to_P3(q), to_P3(r)));
  }

private:
  const Base_functor_ _f;
  const Point_3& _center;
  const FT _radius;
};

// Adaptor for calling the functors with the points projected on the sphere
template <typename LK, typename Functor_>
class Functor_geo_projection_adaptor
  : public Functor_
{
  typedef Functor_                                                   Functor;
  typedef Functor_                                                   Base;

  typedef typename LK::FT                                            FT;
  typedef typename LK::Point_3                                       Point_3;
  typedef Geographical_coordinates<LK>                               Point;

  typedef Construct_Cartesian_coordinates<LK>                        To_Point_3;

public:
  Functor_geo_projection_adaptor(const Functor& f, const Point_3& center, const FT radius)
    : Base(f), to_p3(center, radius)
  { }

public:
  using Base::operator();

  decltype(auto) operator()(const Point& p0, const Point& p1)
  { return Base::operator()(to_p3(p0), to_p3(p1)); }

  decltype(auto) operator()(const Point& p0, const Point& p1, const Point& p2)
  { return Base::operator()(to_p3(p0), to_p3(p1), to_p3(p2)); }

  decltype(auto) operator ()(const Point& p0, const Point& p1, const Point& p2, const Point& p3)
  { return Base::operator()(to_p3(p0), to_p3(p1), to_p3(p2), to_p3(p3)); }

  decltype(auto) operator()(const Point& p0, const Point& p1, const Point& p2, const Point& p3, const Point& p4)
  { return Base::operator()(to_p3(p0), to_p3(p1), to_p3(p2), to_p3(p3), to_p3(p4)); }

private:
  const To_Point_3 to_p3;
};

} // namespace internal

template <typename LK>
class Geographical_coordinates_traits_2
  : public Delaunay_triangulation_on_sphere_traits_2<LK>
{
protected:
  typedef Delaunay_triangulation_on_sphere_traits_2<LK>                        Base;
  typedef Geographical_coordinates_traits_2<LK>                                Self;

public:
  typedef typename LK::FT                                                      FT;
  typedef typename LK::Point_3                                                 Point_3;
  typedef Geographical_coordinates<LK>                                         Point_on_sphere_2;

  // predicates
  typedef internal::Functor_geo_projection_adaptor<LK, typename Base::Collinear_are_strictly_ordered_on_great_circle_2> Collinear_are_strictly_ordered_on_great_circle_2;
  typedef internal::Functor_geo_projection_adaptor<LK, typename Base::Compare_on_sphere_2> Compare_on_sphere_2;
  typedef internal::Functor_geo_projection_adaptor<LK, typename Base::Equal_on_sphere_2> Equal_on_sphere_2;
  typedef internal::Functor_geo_projection_adaptor<LK, typename Base::Orientation_on_sphere_2> Orientation_on_sphere_2;
  typedef internal::Functor_geo_projection_adaptor<LK, typename Base::Side_of_oriented_circle_on_sphere_2> Side_of_oriented_circle_on_sphere_2;

  // constructions
  typedef internal::Construct_geographical_coordinates<LK>                     Construct_point_on_sphere_2;
  typedef internal::Construct_Cartesian_coordinates<LK>                        Construct_point_3;
  typedef internal::Functor_geo_projection_adaptor<LK, typename Base::Construct_segment_3> Construct_segment_3;
  typedef internal::Functor_geo_projection_adaptor<LK, typename Base::Construct_triangle_3> Construct_triangle_3;

  typedef internal::Functor_geo_projection_adaptor<LK, typename Base::Construct_arc_on_sphere_2> Construct_arc_on_sphere_2;
  typedef internal::Construct_circumcenter_on_geo_sphere_2<LK, typename Base::Construct_circumcenter_on_sphere_2> Construct_circumcenter_on_sphere_2;

public:
  Geographical_coordinates_traits_2(const Point_3& sphere = CGAL::ORIGIN,
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
  { return Construct_point_on_sphere_2(Base::center(), Base::radius()); }

  Construct_point_3
  construct_point_3_object() const
  { return Construct_point_3(Base::center(), Base::radius()); }

  Construct_arc_on_sphere_2
  construct_arc_on_sphere_2_object() const
  { return Construct_arc_on_sphere_2(Base::construct_arc_on_sphere_2_object(),
                                     Base::center(), Base::radius()); }

  Construct_circumcenter_on_sphere_2
  construct_circumcenter_on_sphere_2_object() const
  { return Construct_circumcenter_on_sphere_2(Base::construct_circumcenter_on_sphere_2_object(),
                                              Base::center(), Base::radius()); }

public:
  bool is_on_sphere(const Point_on_sphere_2& /*p*/) const
  {
    return true;
  }

  bool are_points_too_close(const Point_on_sphere_2& p, const Point_on_sphere_2& q) const
  {
    return Base::are_points_too_close(construct_point_3_object()(p), construct_point_3_object()(q));
  }
};

} // namespace CGAL

#endif // CGAL_TRIANGULATION_ON_SPHERE_GEOGRAPHICAL_COORDINATES_TRAITS_2_H
