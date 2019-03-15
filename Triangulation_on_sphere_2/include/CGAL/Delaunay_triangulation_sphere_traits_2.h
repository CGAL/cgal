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

template < typename K >
class Power_test_2
{
public:
  typedef typename K::Point_2               Point_2;
  typedef typename K::Oriented_side         Oriented_side;
  typedef typename K::Comparison_result     Comparison_result;

  Power_test_2(const Point_2& sphere);

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

public:
  typedef Oriented_side result_type;

protected:
  Point_2 _sphere;
};

template < typename K >
Power_test_2<K>::
Power_test_2(const Point_2& sphere)
  : _sphere(sphere)
{ }

template < typename K >
class Orientation_sphere_1
{
public:
  typedef typename K::Point_2                  Point_2;
  typedef typename K::Comparison_result        Comparison_result;
  typedef Comparison_result                    result_type;

  Orientation_sphere_1(const Point_2& sphere);

  Comparison_result operator()(const Point_2& p, const Point_2& q) const
  {
    return coplanar_orientation(_sphere, p, q);
  }

  Comparison_result operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
  {
    return coplanar_orientation(p, q, r, _sphere);
  }

  Comparison_result operator()(const Point_2& p, const Point_2& q, const Point_2& r, const Point_2& s) const
  {
    return coplanar_orientation(p, q, r, s);
  }

protected :
  Point_2 _sphere;
};

template < typename K >
Orientation_sphere_1<K>::
Orientation_sphere_1(const Point_2& sphere)
  : _sphere(sphere)
{ }

template < typename K >
class Orientation_sphere_2
{
public:
  typedef typename K::Point_2                  Point_2;
  typedef typename K::Comparison_result        Comparison_result;

  typedef Comparison_result                    result_type;

  Orientation_sphere_2(const Point_2& sphere);

  Comparison_result operator()(const Point_2& p, const Point_2& q,
                               const Point_2& r) const
  { return orientation(_sphere, p, q, r); }

  Comparison_result operator()(const Point_2& p, const Point_2& q,
                               const Point_2& r, const Point_2& s) const
  { return orientation(p, q, r, s); }

protected :
  Point_2 _sphere;
};

template < typename K >
Orientation_sphere_2<K>::
Orientation_sphere_2(const Point_2& sphere)
  : _sphere(sphere)
{ }

template < typename K >
class Coradial_sphere_2
{
public:
  typedef typename K::Point_2                  Point_2;
  typedef bool                                 result_type;

  Coradial_sphere_2(const Point_2& sphere);

  bool operator()(const Point_2& p, const Point_2 q) const
  {
    return collinear(_sphere, p, q) &&
           (are_ordered_along_line(_sphere, p, q) ||
            are_ordered_along_line(_sphere, q, p));
  }

protected :
  Point_2 _sphere;
};

template < typename K >
Coradial_sphere_2<K>::
Coradial_sphere_2(const Point_2& sphere)
  : _sphere(sphere)
{ }

template < typename K >
class Inside_cone_2
{
public:
  typedef typename K::Point_2                  Point_2;
  typedef bool                                 result_type;

  Inside_cone_2(const Point_2& sphere);

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

protected :
  Point_2 _sphere;
};

template < typename K >
Inside_cone_2<K>::
Inside_cone_2(const Point_2& sphere)
  : _sphere(sphere)
{ }

template < class R >
class Delaunay_triangulation_sphere_traits_2
  : public R
{
public:
  typedef typename R::Point_3                       Point_2;
  typedef typename R::Point_3                       Weighted_point_2;
  typedef typename R::Ray_3                         Ray_2;
  typedef typename R::Line_3                        Line_2;
  typedef typename R::Construct_ray_3               Construct_ray_3;
  typedef typename R::Construct_circumcenter_3      Construct_circumcenter_3;
  typedef typename R::Construct_bisector_3          Construct_bisector_3;
  typedef typename R::Construct_segment_3           Construct_segment_3;

  typedef Delaunay_triangulation_sphere_traits_2<R> Self;
  typedef CGAL::Power_test_2<Self>                  Power_test_2;
  typedef CGAL::Orientation_sphere_2<Self>          Orientation_2;
  typedef CGAL::Coradial_sphere_2<Self>             Coradial_sphere_2;
  typedef CGAL::Inside_cone_2<Self>                 Inside_cone_2;
  //typedef CGAL::Orientation_sphere_1<Self>          Orientation_1;
  typedef typename R::Coplanar_orientation_3        Orientation_1;
  typedef typename R::Compute_squared_distance_3    Compute_squared_distance_2;
  typedef typename R::Compare_xyz_3                 Compare_xyz_3;

  typedef boost::true_type                          requires_test;

  double _radius;

  void set_radius(double a) { _radius = a; }

  Delaunay_triangulation_sphere_traits_2(const Point_2& sphere = Point_2(0,0,0));
  Compare_xyz_3
  compare_xyz_3_object() const
  { return Compare_xyz_3(); }

  Compute_squared_distance_2
  compute_squared_distance_2_object() const
  { return Compute_squared_distance_2(); }

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

  Construct_ray_3
  construct_ray_2_object() const
  { return Construct_ray_3(); }

  Construct_circumcenter_3
  construct_circumcenter_2_object() const
  { return Construct_circumcenter_3(); }

  Construct_segment_3
  construct_segment_2_object() const
  { return Construct_segment_3(); }

  Compute_squared_distance_2
  compute_squared_distance_3_object() const
  { return Compute_squared_distance_2(); }

protected :
  Point_2 _sphere;
};

template < class R >
Delaunay_triangulation_sphere_traits_2<R> ::
Delaunay_triangulation_sphere_traits_2(const Point_2& sphere)
  : _sphere(sphere)
{ }

} // namespace CGAL

#endif // CGAL_Reg_TRIANGULATION_SPHERE_TRAITS_2_H
