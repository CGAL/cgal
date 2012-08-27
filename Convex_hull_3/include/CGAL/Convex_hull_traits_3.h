// Copyright (c) 2001,2011  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Susan Hert <hert@mpi-sb.mpg.de>

#ifndef CGAL_CONVEX_HULL_TRAITS_3_H
#define CGAL_CONVEX_HULL_TRAITS_3_H

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Convex_hull_face_base_2.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Projection_traits_xz_3.h>
#include <CGAL/Projection_traits_yz_3.h>
#include <list>

namespace CGAL {
template < class R_ >
class Point_triple
{
protected:
  typedef typename R_::FT                   FT;
  typedef typename R_::Point_3              Point_3;
  typedef typename R_::Vector_3             Vector_3;

  Point_3  p_,  q_,  r_;
public:
  typedef R_                                     R;

  Point_triple() {}

  Point_triple(const Point_3 &p, const Point_3 &q, const Point_3 &r)
    : p_(p), q_(q), r_(r) 
  {}

  const Point_3& p() const { return p_; }
  const Point_3& q() const { return q_; }
  const Point_3& r() const { return r_; }

};

template <class K>
class Point_triple_has_on_positive_side_3 {

public:
    typedef typename K::Point_3 Point_3;
    typedef typename K::Plane_3 Plane_3;
  bool
    operator()( const Plane_3& pl, const Point_3& p) const
    {
      typename K::Orientation_3 o; 
      return ( o(pl.p(), pl.q(), pl.r(), p) == CGAL::POSITIVE );
    }
};
template <class K, class OldK>
class Point_triple_construct_orthogonal_vector_3 
{
public:
  
    typedef typename K::Vector_3 Vector_3;
    typedef typename K::Plane_3 Plane_3;
  
  Vector_3 operator()(const Plane_3& plane) const
  {
    typename OldK::Construct_orthogonal_vector_3 
      construct_orthogonal_vector_3;
    return construct_orthogonal_vector_3(plane.p(), plane.q(), plane.r());
  }
};


template <class K>
class Point_triple_oriented_side_3 
{
public:
  
    typedef typename K::Point_3 Point_3;
    typedef typename K::Plane_3 Plane_3;
    typedef Oriented_side    result_type;

  result_type
    operator()( const Plane_3& pl, const Point_3& p) const
    {
      typename K::Orientation_3 o; 
      Orientation ori = o(pl.p(), pl.q(), pl.r(), p) ;
      if(ori > 0) return ON_POSITIVE_SIDE;
      if(ori < 0) return ON_NEGATIVE_SIDE;
      return ON_ORIENTED_BOUNDARY;
    }
};

template <typename K, typename OldK>
class Point_triple_less_signed_distance_to_plane_3 
{  
public:
    typedef typename K::Point_3 Point_3;
    typedef typename K::Plane_3 Plane_3;

    typedef bool             result_type;

    bool
    operator()( const Plane_3& h, const Point_3& p, const Point_3& q) const
    {
      const Point_3& hp = h.p();
      const Point_3& hq = h.q();
      const Point_3& hr = h.r();
      //typename OldK::Less_signed_distance_to_plane_3 
      //	less_signed_distance_to_plane_3;
      // return less_signed_distance_to_plane_3(hp, hq, hr, p, q);
      return has_smaller_signed_dist_to_planeC3(hp.x(), hp.y(), hp.z(),
                                                hq.x(), hq.y(), hq.z(),
                                                hr.x(), hr.y(), hr.z(),
                                                p.x(), p.y(), p.z(),
                                                q.x(), q.y(), q.z());
      
    }
  };



template <class T>
class Max_coordinate_3 
{
public:

    int operator()(const T& v)
    {
      if (CGAL_NTS abs(v.x()) >= CGAL_NTS abs(v.y()))
      {
         if (CGAL_NTS abs(v.x()) >= CGAL_NTS abs(v.z())) return 0;
         return 2;
      }
      else
      {
         if (CGAL_NTS abs(v.y()) >= CGAL_NTS abs(v.z())) return 1;
         return 2;
      }
    }
};


template <typename GT>
struct GT3_for_CH3 {
  typedef typename GT::Point_3 Point_2;
};



template <class R_>
class Convex_hull_traits_3 
{
 public:  
  typedef R_                                     R;
  typedef Convex_hull_traits_3<R>                Self;
  typedef typename R::Point_3                    Point_3;
  typedef typename R::Segment_3                  Segment_3;
  typedef typename R::Triangle_3                 Triangle_3;
  typedef Point_triple<R>                        Plane_3;
  typedef typename R::Vector_3                   Vector_3;

  typedef CGAL::Polyhedron_3<R>                  Polyhedron_3;

  typedef typename R::Construct_segment_3        Construct_segment_3;
  typedef typename R::Construct_ray_3            Construct_ray_3;

  class Construct_plane_3 {
  public:
    Plane_3 operator ()(const Point_3& p, const Point_3& q, const Point_3& r)
    {
      return Plane_3(p,q,r);
    }
  };

  typedef typename R::Construct_vector_3         Construct_vector_3;
  typedef typename R::Construct_triangle_3       Construct_triangle_3;
  typedef typename R::Construct_centroid_3       Construct_centroid_3;
  typedef Point_triple_construct_orthogonal_vector_3<Self, R>
                                                 Construct_orthogonal_vector_3;

  typedef typename R::Equal_3                    Equal_3;
  typedef typename R::Orientation_3              Orientation_3;
  typedef typename R::Collinear_3                Collinear_3;
  typedef typename R::Coplanar_3                 Coplanar_3;
  typedef typename R::Less_distance_to_point_3   Less_distance_to_point_3;
  typedef Point_triple_has_on_positive_side_3<Self>     Has_on_positive_side_3;

  typedef  Point_triple_less_signed_distance_to_plane_3<Self, R>
                                               Less_signed_distance_to_plane_3;

  // required for degenerate case of all points coplanar
  typedef CGAL::Max_coordinate_3<Vector_3>             Max_coordinate_3;

  // for postcondition checking 
  typedef typename R::Ray_3                      Ray_3; 

  typedef typename R::Has_on_3                   Has_on_3;
  typedef Point_triple_oriented_side_3<Self>     Oriented_side_3;
  typedef typename R::Intersect_3                Intersect_3;
  typedef typename R::Do_intersect_3             Do_intersect_3;

  Construct_segment_3
  construct_segment_3_object() const
  { return Construct_segment_3(); }

  Construct_ray_3
  construct_ray_3_object() const
  { return Construct_ray_3(); }

  Construct_plane_3
  construct_plane_3_object() const
  { return Construct_plane_3(); }

  Construct_triangle_3
  construct_triangle_3_object() const
  { return Construct_triangle_3(); }

  Construct_vector_3
  construct_vector_3_object() const
  { return Construct_vector_3(); }

  Construct_centroid_3  
  construct_centroid_3_object() const
  { return Construct_centroid_3(); }

  Construct_orthogonal_vector_3
  construct_orthogonal_vector_3_object() const
  { return Construct_orthogonal_vector_3(); }

  Collinear_3
  collinear_3_object() const
  { return Collinear_3(); }

  Coplanar_3
  coplanar_3_object() const
  { return Coplanar_3(); }
 
  Has_on_3
  has_on_3_object() const
  { return Has_on_3(); }

  Less_distance_to_point_3
  less_distance_to_point_3_object() const
  { return Less_distance_to_point_3(); }

  Has_on_positive_side_3
  has_on_positive_side_3_object() const
  { return Has_on_positive_side_3(); }

  Oriented_side_3
  oriented_side_3_object() const
  { return Oriented_side_3(); }

  Equal_3
  equal_3_object() const
  { return Equal_3(); }

  Intersect_3
  intersect_3_object() const
  { return Intersect_3(); }

  Do_intersect_3
  do_intersect_3_object() const
  { return Do_intersect_3(); }

  Less_signed_distance_to_plane_3  
  less_signed_distance_to_plane_3_object() const
  { return Less_signed_distance_to_plane_3(); }

  Max_coordinate_3  
  max_coordinate_3_object() const
  { return Max_coordinate_3(); }
};

} // namespace CGAL

#endif // CGAL_CONVEX_HULL_TRAITS_3_H
