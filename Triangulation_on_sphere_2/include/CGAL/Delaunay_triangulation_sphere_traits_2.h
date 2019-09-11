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

#ifndef CGAL_DELAUNAY_TRIANGULATION_SPHERE_TRAITS_2_H
#define CGAL_DELAUNAY_TRIANGULATION_SPHERE_TRAITS_2_H

#include <CGAL/enum.h>

namespace CGAL {

template <typename TTraits>
class Power_test_2
{
public:
  typedef typename TTraits::Point_2               Point_2;
  typedef typename TTraits::Oriented_side         Oriented_side;
  typedef typename TTraits::Comparison_result     Comparison_result;
  typedef Oriented_side                           result_type;

  Power_test_2(const Point_2& sphere) : _sphere(sphere) { }

  Oriented_side operator()(const Point_2& p,
                           const Point_2& q,
                           const Point_2& r,
                           const Point_2& s) const
  {
    return orientation(p, q, r, s);
  }

  Oriented_side operator()(const Point_2& p,
                           const Point_2& q,
                           const Point_2& r) const
  {
    return -coplanar_orientation(p, q,_sphere, r);
  }

  Oriented_side operator()(const Point_2& p,
                           const Point_2& q) const
  {
    Comparison_result pq = compare_xyz(p, q);

    if(pq == EQUAL)
      return ON_ORIENTED_BOUNDARY;

    Comparison_result sq = compare_xyz(_sphere, q);
    if(pq == sq)
      return ON_POSITIVE_SIDE;

    return ON_NEGATIVE_SIDE;
  }

protected:
  const Point_2& _sphere;
};

template <typename TTraits>
class Orientation_sphere_2
{
public:
  typedef typename TTraits::Point_2                  Point_2;
  typedef typename TTraits::Comparison_result        Comparison_result;
  typedef Comparison_result                          result_type;

  Orientation_sphere_2(const Point_2& sphere) : _sphere(sphere) { }

  Comparison_result operator()(const Point_2& p, const Point_2& q,
                               const Point_2& r) const
  { return orientation(_sphere, p, q, r); }

  Comparison_result operator()(const Point_2& p, const Point_2& q,
                               const Point_2& r, const Point_2& s) const
  { return orientation(p, q, r, s); }

protected:
  const Point_2& _sphere;
};

template <typename TTraits>
class Coradial_sphere_2
{
public:
  typedef typename TTraits::Point_2                  Point_2;
  typedef bool                                       result_type;

  Coradial_sphere_2(const Point_2& sphere) : _sphere(sphere) { }

  bool operator()(const Point_2& p, const Point_2 q) const
  {
    return collinear(_sphere, p, q) &&
           (are_ordered_along_line(_sphere, p, q) ||
            are_ordered_along_line(_sphere, q, p));
  }

protected:
  const Point_2& _sphere;
};

template <typename TTraits>
class Inside_cone_2
{
public:
  typedef typename TTraits::Point_2                  Point_2;
  typedef bool                                       result_type;

  Inside_cone_2(const Point_2& sphere) : _sphere(sphere) { }

  bool operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
  {
    if(collinear(_sphere, p, r) ||
       collinear(_sphere, q, r) ||
       orientation(_sphere, p, q, r) != COLLINEAR)
      return false;

    if(collinear(_sphere, p, q))
      return true;

    return coplanar_orientation(_sphere, p, q, r) ==
             (POSITIVE == coplanar_orientation(_sphere, q, p, r));
  }

protected:
  const Point_2& _sphere;
};

template <typename K>
class Delaunay_triangulation_sphere_traits_2
  : public K
{
  typedef K                                         Base;
  typedef Delaunay_triangulation_sphere_traits_2<K> Self;

public:
  typedef typename K::FT                            FT;
  typedef typename K::Point_3                       Point_2;
  typedef typename K::Point_3                       Point;
  typedef typename K::Point_3                       Weighted_point_2; // @fixme
  typedef typename K::Ray_3                         Ray_2;
  typedef typename K::Line_3                        Line_2;

  typedef typename K::Compare_xyz_3                 Compare_xyz_3;
  typedef typename K::Construct_bisector_3          Construct_bisector_3;
  typedef typename K::Construct_circumcenter_3      Construct_circumcenter_3;
  typedef typename K::Construct_ray_3               Construct_ray_3;
  typedef typename K::Construct_segment_3           Construct_segment_3;

  typedef typename K::Compute_squared_distance_3    Compute_squared_distance_2;

  typedef CGAL::Coradial_sphere_2<Self>             Coradial_sphere_2;
  typedef CGAL::Inside_cone_2<Self>                 Inside_cone_2;
  typedef CGAL::Orientation_sphere_2<Self>          Orientation_2;
  typedef CGAL::Power_test_2<Self>                  Power_test_2;

  typedef typename K::Coplanar_orientation_3        Orientation_1;

  Delaunay_triangulation_sphere_traits_2(const Point_2& sphere = CGAL::ORIGIN,
                                         const FT radius = 1,
                                         const K& k = K())
    : Base(k), _sphere(sphere), _radius(radius)
  {
    initialize_bounds();
  }

private:
  void initialize_bounds()
  {
    const FT minDist = _radius * std::pow(2, -23); // @fixme
    const FT minRadius = _radius * (1 - std::pow(2, -50));
    const FT maxRadius = _radius * (1 + std::pow(2, -50));
    _minDistSquared = CGAL::square(minDist);
    _minRadiusSquared = CGAL::square(minRadius);
    _maxRadiusSquared = CGAL::square(maxRadius);
  }

public:
  Point_2 center() const { return _sphere; }
  void set_center(const Point_2& sphere) { _sphere = sphere; }
  void set_radius(const FT radius) { _radius = radius; initialize_bounds(); }

  bool is_on_sphere(const Point_2& p) const
  {
    const FT sq_dist = Base::compute_squared_distance_3_object()(p, _sphere);
    return (_minRadiusSquared < sq_dist && sq_dist < _maxRadiusSquared);
  }

  bool are_points_too_close(const Point_2& p, const Point_2& q) const
  {
    return (Base::compute_squared_distance_3_object()(p, q) <= _minDistSquared);
  }

  Compute_squared_distance_2
  compute_squared_distance_3_object() const
  { return Base::compute_squared_distance_3_object(); }

  Construct_ray_3
  construct_ray_2_object() const
  { return Base::construct_ray_3_object(); }

  Construct_circumcenter_3
  construct_circumcenter_2_object() const
  { return Base::construct_circumcenter_3_object(); }

  Construct_segment_3
  construct_segment_2_object() const
  { return Base::construct_segment_3_object(); }

  Orientation_2
  orientation_2_object() const
  { return Orientation_2(_sphere); }

  Orientation_1
  orientation_1_object() const
  { return Orientation_1(); }

  Power_test_2
  power_test_2_object() const
  { return Power_test_2(_sphere); }

  Coradial_sphere_2
  coradial_sphere_2_object() const
  { return Coradial_sphere_2(_sphere); }

  Inside_cone_2
  inside_cone_2_object() const
  { return Inside_cone_2(_sphere); }

protected:
  Point_2 _sphere;
  FT _radius;

  FT _minDistSquared; // minimal distance of two points to each other
  FT _minRadiusSquared; // minimal distance of a point from center of the sphere
  FT _maxRadiusSquared; // maximal distance of a point from center of the sphere
};

} // namespace CGAL

#endif // CGAL_DELAUNAY_TRIANGULATION_SPHERE_TRAITS_2_H
