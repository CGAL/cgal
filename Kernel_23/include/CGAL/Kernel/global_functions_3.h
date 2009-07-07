// Copyright (c) 2003-2004  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
// Author(s)     : Sylvain Pion
 
#ifndef CGAL_KERNEL_GLOBAL_FUNCTIONS_3_H
#define CGAL_KERNEL_GLOBAL_FUNCTIONS_3_H

#include <CGAL/user_classes.h>
#include <CGAL/Kernel/global_functions_internal_3.h>
#include <CGAL/Kernel/mpl.h>

// Generic functions calling the kernel functor.
// See comments in CGAL/Kernel/global_functions_2.h.

CGAL_BEGIN_NAMESPACE

template <typename K>
inline
Angle
angle(const Point_3<K> &p, const Point_3<K> &q, const Point_3<K> &r)
{
  return CGALi::angle(p, q, r, K());
}

template < typename K >
inline
typename K::Boolean
are_negative_oriented(const Point_3<K>& p, const Point_3<K>& q,
                      const Point_3<K>& r, const Point_3<K>& s)
{
  return CGALi::are_negative_oriented(p, q, r, s, K());
}

template < class K >
inline
typename K::Boolean
are_ordered_along_line(const Point_3<K> &p,
                       const Point_3<K> &q,
                       const Point_3<K> &r)
{
  return CGALi::are_ordered_along_line(p, q, r, K());
}

template < typename K >
inline
typename K::Boolean
are_positive_oriented(const Point_3<K>& p, const Point_3<K>& q,
                      const Point_3<K>& r, const Point_3<K>& s)
{
  return CGALi::are_positive_oriented(p, q, r, s, K());
}

template < class K >
inline
typename K::Boolean
are_strictly_ordered_along_line(const Point_3<K> &p,
                                const Point_3<K> &q,
                                const Point_3<K> &r)
{
  return CGALi::are_strictly_ordered_along_line(p, q, r, K());
}

template < class K >
inline
typename K::Point_3
barycenter(const Point_3<K> &p1, const typename K::FT& w1,
           const Point_3<K> &p2)
{
  return CGALi::barycenter(p1, w1, p2, K());
}

template < class K >
inline
typename K::Point_3
barycenter(const Point_3<K> &p1, const typename K::FT& w1,
           const Point_3<K> &p2, const typename K::FT& w2)
{
  return CGALi::barycenter(p1, w1, p2, w2, K());
}

template < class K >
inline
typename K::Point_3
barycenter(const Point_3<K> &p1, const typename K::FT& w1,
           const Point_3<K> &p2, const typename K::FT& w2,
           const Point_3<K> &p3)
{
  return CGALi::barycenter(p1, w1, p2, w2, p3, K());
}

template < class K >
inline
typename K::Point_3
barycenter(const Point_3<K> &p1, const typename K::FT& w1,
           const Point_3<K> &p2, const typename K::FT& w2,
           const Point_3<K> &p3, const typename K::FT& w3)
{
  return CGALi::barycenter(p1, w1, p2, w2, p3, w3, K());
}

template < class K >
inline
typename K::Point_3
barycenter(const Point_3<K> &p1, const typename K::FT& w1,
           const Point_3<K> &p2, const typename K::FT& w2,
           const Point_3<K> &p3, const typename K::FT& w3,
           const Point_3<K> &p4)
{
  return CGALi::barycenter(p1, w1, p2, w2, p3, w3, p4, K());
}

template < class K >
inline
typename K::Point_3
barycenter(const Point_3<K> &p1, const typename K::FT& w1,
           const Point_3<K> &p2, const typename K::FT& w2,
           const Point_3<K> &p3, const typename K::FT& w3,
           const Point_3<K> &p4, const typename K::FT& w4)
{
  return CGALi::barycenter(p1, w1, p2, w2, p3, w3, p4, w4, K());
}

template <typename K>
inline
typename K::Plane_3
bisector(const Point_3<K> &p, const Point_3<K> &q)
{
  return CGALi::bisector(p, q, K());
}

template <typename K>
inline
typename K::Plane_3
bisector(const Plane_3<K> &h1, const Plane_3<K> &h2)
{
  return CGALi::bisector(h1, h2, K());
}

template < class K >
inline
Point_3<K>
centroid(const Point_3<K> &p, const Point_3<K> &q,
         const Point_3<K> &r, const Point_3<K> &s)
{
  return CGALi::centroid(p, q, r, s, K());
}

template < class K >
inline
Point_3<K>
centroid(const Point_3<K> &p, const Point_3<K> &q, const Point_3<K> &r)
{
  return CGALi::centroid(p, q, r, K());
}

template < class K >
inline
Point_3<K>
centroid(const Tetrahedron_3<K> &t)
{
  return CGALi::centroid(t, K());
}

template < class K >
inline
Point_3<K>
centroid(const Triangle_3<K> &t)
{
  return CGALi::centroid(t, K());
}

template < class K >
inline
typename K::Point_3
circumcenter(const Point_3<K> &p,
             const Point_3<K> &q)
{
  return CGALi::circumcenter(p, q, K());
}

template < class K >
inline
typename K::Point_3
circumcenter(const Point_3<K> &p,
             const Point_3<K> &q,
             const Point_3<K> &r)
{
  return CGALi::circumcenter(p, q, r, K());
}

template < class K >
inline
typename K::Point_3
circumcenter(const Point_3<K> &p, const Point_3<K> &q,
             const Point_3<K> &r, const Point_3<K> &s)
{
  return CGALi::circumcenter(p, q, r, s, K());
}

template < class K >
inline
typename K::Point_3
circumcenter(const Tetrahedron_3<K> &t)
{
  return CGALi::circumcenter(t, K());
}

template < class K >
inline
typename K::Point_3
circumcenter(const Triangle_3<K> &t)
{
  return CGALi::circumcenter(t, K());
}

template < class K >
inline
typename K::Boolean
collinear(const Point_3<K> &p, const Point_3<K> &q, const Point_3<K> &r)
{
  return CGALi::collinear(p, q, r, K());
}

template < class K >
inline
typename K::Boolean
collinear_are_ordered_along_line(const Point_3<K> &p,
                                 const Point_3<K> &q,
                                 const Point_3<K> &r)
{
  return CGALi::collinear_are_ordered_along_line(p, q, r, K());
}

template < class K >
inline
typename K::Boolean
collinear_are_strictly_ordered_along_line(const Point_3<K> &p,
                                          const Point_3<K> &q,
                                          const Point_3<K> &r)
{
  return CGALi::collinear_are_strictly_ordered_along_line(p, q, r, K());
}

template < class K >
inline
typename K::Comparison_result
compare_distance_to_point(const Point_3<K> &p,
                          const Point_3<K> &q,
                          const Point_3<K> &r)
{
  return CGALi::compare_distance_to_point(p, q, r, K());
}

template < class K >
inline
typename K::Comparison_result
compare_squared_distance(const Point_3<K> &p,
                         const Point_3<K> &q,
                         const typename K::FT &d2)
{
  return CGALi::compare_squared_distance(p, q, d2, K());
}

template < class K >
inline
typename K::Comparison_result
compare_squared_radius(const Point_3<K> &p,
		       const Point_3<K> &q,
		       const typename K::FT &sr)
{
  return CGALi::compare_squared_radius(p, q, sr, K());
}

template < class K >
inline
typename K::Comparison_result
compare_squared_radius(const Point_3<K> &p,
		       const Point_3<K> &q,
		       const Point_3<K> &r,
		       const typename K::FT &sr)
{
  return CGALi::compare_squared_radius(p, q, r, sr, K());
}

template < class K >
inline
typename K::Comparison_result
compare_squared_radius(const Point_3<K> &p,
		       const Point_3<K> &q,
		       const Point_3<K> &r,
		       const Point_3<K> &s,
		       const typename K::FT &sr)
{
  return CGALi::compare_squared_radius(p, q, r, s, sr, K());
}

template < class K >
inline
typename K::Comparison_result
compare_lexicographically_xyz(const Point_3<K> &p,
                              const Point_3<K> &q)
{
  return CGALi::compare_lexicographically_xyz(p, q, K());
}

template < class K >
inline
typename K::Comparison_result
compare_signed_distance_to_plane(const Plane_3<K> &h,
				 const Point_3<K> &p,
				 const Point_3<K> &q)
{ 
  return CGALi::compare_signed_distance_to_plane(h, p, q, K());
}

template < class K >
inline
typename K::Comparison_result
compare_signed_distance_to_plane(const Point_3<K> &hp,
				 const Point_3<K> &hq,
				 const Point_3<K> &hr,
				 const Point_3<K> &p,
				 const Point_3<K> &q)
{ 
  return CGALi::compare_signed_distance_to_plane(hp, hq, hr, p, q, K());
}

template < class K >
inline
typename K::Comparison_result
compare_x(const Point_3<K> &p, const Point_3<K> &q)
{ 
  return CGALi::compare_x(p, q, K());
}

template < class K >
inline
typename K::Comparison_result
compare_y(const Point_3<K> &p, const Point_3<K> &q)
{ 
  return CGALi::compare_y(p, q, K());
}

template < class K >
inline
typename K::Comparison_result
compare_z(const Point_3<K> &p, const Point_3<K> &q)
{ 
  return CGALi::compare_z(p, q, K());
}

template < class K >
inline
typename K::Comparison_result
compare_xyz(const Point_3<K> &p, const Point_3<K> &q)
{ 
  return CGALi::compare_xyz(p, q, K());
}

template < class K >
inline
typename K::Boolean
coplanar(const Point_3<K> &p, const Point_3<K> &q,
         const Point_3<K> &r, const Point_3<K> &s)
{
  return CGALi::coplanar(p, q, r, s, K());
}


template < class K >
inline
typename K::Orientation
coplanar_orientation(const Point_3<K> &p,
                     const Point_3<K> &q,
                     const Point_3<K> &r,
                     const Point_3<K> &s)
{
  return CGALi::coplanar_orientation(p, q, r, s, K());
}

template < class K >
inline
typename K::Orientation
coplanar_orientation(const Point_3<K> &p,
                     const Point_3<K> &q,
                     const Point_3<K> &r)
{
  return CGALi::coplanar_orientation(p, q, r, K());
}

template < class K >
inline
typename K::Bounded_side
coplanar_side_of_bounded_circle(const Point_3<K> &p,
                                const Point_3<K> &q,
                                const Point_3<K> &r,
                                const Point_3<K> &t)
{
  return CGALi::coplanar_side_of_bounded_circle(p, q, r, t, K());
}

template < class K >
inline
typename K::Vector_3
cross_product(const Vector_3<K> &v, const Vector_3<K> &w)
{
  return CGALi::cross_product(v, w, K());
}

template < class K >
inline
typename K::FT
determinant(const Vector_3<K> &v0, const Vector_3<K> &v1,
            const Vector_3<K> &v2)
{
  return CGALi::determinant(v0, v1, v2, K());
}

template < class K >
inline
typename K::Line_3
equidistant_line(const Point_3<K> &p, const Point_3<K> &q, const Point_3<K> &r)
{
  return CGALi::equidistant_line(p, q, r, K());
}

template < class K >
inline
typename K::Boolean
has_larger_distance_to_point(const Point_3<K> &p,
			     const Point_3<K> &q,
			     const Point_3<K> &r)
{
  return CGALi::has_larger_distance_to_point(p, q, r, K());
}

template < class K >
inline
typename K::Boolean
has_larger_signed_distance_to_plane(const Point_3<K> &hp,
				    const Point_3<K> &hq,
				    const Point_3<K> &hr,
				    const Point_3<K> &p,
				    const Point_3<K> &q)
{ 
  return CGALi::has_larger_signed_distance_to_plane(hp, hq, hr, p, q, K());
}

template < class K >
inline
typename K::Boolean
has_larger_signed_distance_to_plane(const Plane_3<K> &h,
				    const Point_3<K> &p,
				    const Point_3<K> &q)
{ 
  return CGALi::has_larger_signed_distance_to_plane(h, p, q, K());
}

template < class K >
inline
typename K::Boolean
has_smaller_distance_to_point(const Point_3<K> &p,
                              const Point_3<K> &q,
                              const Point_3<K> &r)
{
  return CGALi::has_smaller_distance_to_point(p, q, r, K());
}

template < class K >
inline
typename K::Boolean
has_smaller_signed_distance_to_plane(const Point_3<K> &hp,
                                     const Point_3<K> &hq,
                                     const Point_3<K> &hr,
                                     const Point_3<K> &p,
                                     const Point_3<K> &q)
{ 
  return CGALi::has_smaller_signed_distance_to_plane(hp, hq, hr, p, q, K());
}

template < class K >
inline
typename K::Boolean
has_smaller_signed_distance_to_plane(const Plane_3<K> &h,
                                     const Point_3<K> &p,
                                     const Point_3<K> &q)
{ 
  return CGALi::has_smaller_signed_distance_to_plane(h, p, q, K());
}

template < class K >
inline
typename K::Boolean
less_x(const Point_3<K> &p, const Point_3<K> &q)
{ 
  return CGALi::less_x(p, q, K());
}

template < class K >
inline
typename K::Boolean
less_y(const Point_3<K> &p, const Point_3<K> &q)
{ 
  return CGALi::less_y(p, q, K());
}

template < class K >
inline
typename K::Boolean
less_z(const Point_3<K> &p, const Point_3<K> &q)
{ 
  return CGALi::less_z(p, q, K());
}

template < class K >
inline
typename K::Boolean
lexicographically_xyz_smaller(const Point_3<K> &p, const Point_3<K> &q)
{
  return CGALi::lexicographically_xyz_smaller(p, q, K());
}

template < class K >
inline
typename K::Boolean
lexicographically_xyz_smaller_or_equal(const Point_3<K> &p,
                                       const Point_3<K> &q)
{
  return CGALi::lexicographically_xyz_smaller_or_equal(p, q, K());
}

template < class K >
inline
typename K::Point_3
midpoint(const Point_3<K> &p, const Point_3<K> &q)
{
  return CGALi::midpoint(p, q, K());
}

template < class K >
inline
typename K::Point_3
max_vertex(const Iso_cuboid_3<K> &ic)
{
  return CGALi::max_vertex(ic, K());
}

template < class K >
inline
typename K::Point_3
min_vertex(const Iso_cuboid_3<K> &ic)
{
  return CGALi::min_vertex(ic, K());
}

template < class K >
inline
typename K::Vector_3
normal(const Point_3<K> &p, const Point_3<K> &q, const Point_3<K> &r)
{
  return CGALi::normal(p, q, r, K());
}

// FIXME TODO : what to do with the operators ?
template < class K >
inline
typename K::Boolean
operator==(const Point_3<K>& p, const Point_3<K>& q)
{ return K().equal_3_object()(p, q); }

template < class K >
inline
typename K::Boolean
operator!=(const Point_3<K>& p, const Point_3<K>& q)
{ return ! (p == q); }

template < class K >
inline
typename K::Boolean
operator==(const Point_3<K>& p, const Origin& o)
{ return K().equal_3_object()(p, Point_3<K>(o)); }

template < class K >
inline
typename K::Boolean
operator!=(const Point_3<K>& p, const Origin& o)
{ return ! (p == o); }

template < class K >
inline
typename K::Boolean
operator==(const Iso_cuboid_3<K>& p, const Iso_cuboid_3<K>& q)
{ return K().equal_3_object()(p, q); }

template < class K >
inline
typename K::Boolean
operator!=(const Iso_cuboid_3<K>& p, const Iso_cuboid_3<K>& q)
{ return ! (p == q); }

template < class K >
inline
typename K::Boolean
operator==(const Plane_3<K>& p, const Plane_3<K>& q)
{ return K().equal_3_object()(p, q); }

template < class K >
inline
typename K::Boolean
operator!=(const Plane_3<K>& p, const Plane_3<K>& q)
{ return ! (p == q); }

template < class K >
inline
typename K::Boolean
operator==(const Segment_3<K>& p, const Segment_3<K>& q)
{ return K().equal_3_object()(p, q); }

template < class K >
inline
typename K::Boolean
operator!=(const Segment_3<K>& p, const Segment_3<K>& q)
{ return ! (p == q); }

template < class K >
inline
typename K::Boolean
operator==(const Line_3<K>& p, const Line_3<K>& q)
{ return K().equal_3_object()(p, q); }

template < class K >
inline
typename K::Boolean
operator!=(const Line_3<K>& p, const Line_3<K>& q)
{ return ! (p == q); }

template < class K >
inline
typename K::Boolean
operator==(const Ray_3<K>& p, const Ray_3<K>& q)
{ return K().equal_3_object()(p, q); }

template < class K >
inline
typename K::Boolean
operator!=(const Ray_3<K>& p, const Ray_3<K>& q)
{ return ! (p == q); }

template < class K >
inline
typename K::Boolean
operator==(const Triangle_3<K>& p, const Triangle_3<K>& q)
{ return K().equal_3_object()(p, q); }

template < class K >
inline
typename K::Boolean
operator!=(const Triangle_3<K>& p, const Triangle_3<K>& q)
{ return ! (p == q); }

template < class K >
inline
typename K::Boolean
operator==(const Tetrahedron_3<K>& p, const Tetrahedron_3<K>& q)
{ return K().equal_3_object()(p, q); }

template < class K >
inline
typename K::Boolean
operator!=(const Tetrahedron_3<K>& p, const Tetrahedron_3<K>& q)
{ return ! (p == q); }

template < class K >
inline
typename K::Boolean
operator==(const Direction_3<K>& p, const Direction_3<K>& q)
{ return K().equal_3_object()(p, q); }

template < class K >
inline
typename K::Boolean
operator!=(const Direction_3<K>& p, const Direction_3<K>& q)
{ return ! (p == q); }

template < class K >
inline
typename K::Boolean
operator==(const Sphere_3<K>& p, const Sphere_3<K>& q)
{ return K().equal_3_object()(p, q); }

template < class K >
inline
typename K::Boolean
operator!=(const Sphere_3<K>& p, const Sphere_3<K>& q)
{ return ! (p == q); }

template < class K >
inline
typename K::Boolean
operator==(const Vector_3<K>& p, const Vector_3<K>& q)
{ return K().equal_3_object()(p, q); }

template < class K >
inline
typename K::Boolean
operator!=(const Vector_3<K>& p, const Vector_3<K>& q)
{ return ! (p == q); }

template < class K >
inline
typename K::Boolean
operator==(const Vector_3<K>& p, const Null_vector& o)
{ return K().equal_3_object()(p, Vector_3<K>(o)); }

template < class K >
inline
typename K::Boolean
operator!=(const Vector_3<K>& p, const Null_vector& o)
{ return ! (p == o); }


template < class K >
inline
typename K::Boolean
operator<(const Point_3<K>& p, const Point_3<K>& q)
{ return K().less_xyz_3_object()(p, q); }

template < class K >
inline
typename K::Boolean
operator>(const Point_3<K>& p, const Point_3<K>& q)
{ return K().less_xyz_3_object()(q, p); }

template < class K >
inline
typename K::Boolean
operator<=(const Point_3<K>& p, const Point_3<K>& q)
{ return ! K().less_xyz_3_object()(q, p); }

template < class K >
inline
typename K::Boolean
operator>=(const Point_3<K>& p, const Point_3<K>& q)
{ return ! K().less_xyz_3_object()(p, q); }

template < class K >
inline
typename K::Vector_3
operator*(const typename K::FT &c, const Vector_3<K> &w)
{
  return K().construct_scaled_vector_3_object()(w, c);
}

template < class K >
inline
typename K::Vector_3
operator*(const Vector_3<K> &w, const typename K::FT &c)
{
  return K().construct_scaled_vector_3_object()(w, c);
}

template < class K >
inline
typename K::Vector_3
operator*(const typename First_if_different<typename K::RT,
                                            typename K::FT>::Type &c,
          const Vector_3<K> &w)
{
  return K().construct_scaled_vector_3_object()(w, c);
}

template < class K >
inline
typename K::Vector_3
operator*(const Vector_3<K> &w,
          const typename First_if_different<typename K::RT,
                                            typename K::FT>::Type &c)
{
  return K().construct_scaled_vector_3_object()(w, c);
}

template < class K >
inline
typename K::FT
operator*(const Vector_3<K> &v, const Vector_3<K> &w)
{
  return K().compute_scalar_product_3_object()(v, w);
}


template < class K >
inline
typename K::Point_3
operator+(const Point_3<K> &p, const Vector_3<K> &v)
{
  return K().construct_translated_point_3_object()(p, v);
}

template < class K >
inline
typename K::Point_3
operator+(const Origin &o, const Vector_3<K> &v)
{
  return K().construct_translated_point_3_object()(o, v);
}

template < class K >
inline
typename K::Point_3
operator-(const Point_3<K> &p, const Vector_3<K> &v)
{
  return K().construct_translated_point_3_object()
               (p, K().construct_opposite_vector_3_object()(v));
}

template < class K >
inline
typename K::Point_3
operator-(const Origin &o, const Vector_3<K> &v)
{
  return K().construct_translated_point_3_object()
               (o, K().construct_opposite_vector_3_object()(v));
}

template < class K >
inline
typename K::Vector_3
operator-(const Point_3<K> &p, const Point_3<K> &q)
{
  return K().construct_vector_3_object()(q, p);
}

template < class K >
inline
typename K::Vector_3
operator-(const Point_3<K> &p, const Origin &o)
{
  return K().construct_vector_3_object()(o, p);
}

template < class K >
inline
typename K::Vector_3
operator-(const Origin &o, const Point_3<K> &q)
{
  return K().construct_vector_3_object()(q, o);
}

template <class K >
inline
typename K::Orientation
orientation(const Point_3<K> &p,
            const Point_3<K> &q,
            const Point_3<K> &r,
            const Point_3<K> &s)
{
  return CGALi::orientation(p, q, r, s, K());
}

template <class K >
inline
typename K::Orientation
orientation(const Vector_3<K> &u, const Vector_3<K> &v, const Vector_3<K> &w)
{
  return CGALi::orientation(u, v, w, K());
}

template <class K >
inline
typename K::Vector_3
orthogonal_vector(const Point_3<K>& p,
		  const Point_3<K>& q,
		  const Point_3<K>& r)
{
  return CGALi::orthogonal_vector(p, q, r, K());
}

template <class K >
inline
typename K::Vector_3
orthogonal_vector(const Plane_3<K>& p)
{
  return CGALi::orthogonal_vector(p, K());
}

// parallel() functions are in Kernel/global_functions.h

template <class K>
inline
typename K::Plane_3
radical_plane(const Sphere_3<K> &s1,
              const Sphere_3<K> &s2)
{
  return K().construct_radical_plane_3_object()(s1,s2);
}

template <class K >
inline
typename K::Bounded_side
side_of_bounded_sphere(const Point_3<K> &p,
                       const Point_3<K> &q,
                       const Point_3<K> &test)
{
  return CGALi::side_of_bounded_sphere(p, q, test, K());
}

template <class K >
inline
typename K::Bounded_side
side_of_bounded_sphere(const Point_3<K> &p,
                       const Point_3<K> &q,
                       const Point_3<K> &r,
                       const Point_3<K> &test)
{
  return CGALi::side_of_bounded_sphere(p, q, r, test, K());
}

template <class K >
inline
typename K::Bounded_side
side_of_bounded_sphere(const Point_3<K> &p,
                       const Point_3<K> &q,
                       const Point_3<K> &r,
                       const Point_3<K> &s,
                       const Point_3<K> &test)
{
  return CGALi::side_of_bounded_sphere(p, q, r, s, test, K());
}

template <class K >
inline
typename K::Oriented_side
side_of_oriented_sphere(const Point_3<K> &p,
                        const Point_3<K> &q,
                        const Point_3<K> &r,
                        const Point_3<K> &s,
                        const Point_3<K> &test)
{
  return CGALi::side_of_oriented_sphere(p, q, r, s, test, K());
}

template <typename K>
inline
typename K::FT
squared_area(const Point_3<K> &p, const Point_3<K> &q, const Point_3<K> &r)
{
  return CGALi::squared_area(p, q, r, K());
}

template < class K >
inline
typename K::FT
squared_radius(const Point_3<K> &p, const Point_3<K> &q,
	       const Point_3<K> &r, const Point_3<K> &s)
{
  return CGALi::squared_radius(p, q, r, s, K());
}

template < class K >
inline
typename K::FT
squared_radius(const Point_3<K> &p, const Point_3<K> &q, const Point_3<K> &r)
{
  return CGALi::squared_radius(p, q, r, K());
}

template < class K >
inline
typename K::FT
squared_radius(const Point_3<K> &p, const Point_3<K> &q)
{
  return CGALi::squared_radius(p, q, K());
}

template < class K >
inline
typename K::FT
squared_radius(const Point_3<K> &p)
{
  return CGALi::squared_radius(p, K());
}

template < class K >
inline
typename K::Vector_3
unit_normal(const Point_3<K> &p, const Point_3<K> &q, const Point_3<K> &r)
{
  return CGALi::unit_normal(p, q, r, K());
}

template < class K >
inline
typename K::FT
volume(const Point_3<K> &p, const Point_3<K> &q,
       const Point_3<K> &r, const Point_3<K> &s)
{
  return CGALi::volume(p, q, r, s, K());
}

template < class K >
inline
typename K::Boolean
x_equal(const Point_3<K> &p, const Point_3<K> &q)
{
  return CGALi::x_equal(p, q, K());
}

template < class K >
inline
typename K::Boolean
y_equal(const Point_3<K> &p, const Point_3<K> &q)
{
  return CGALi::y_equal(p, q, K());
}

template < class K >
inline
typename K::Boolean
z_equal(const Point_3<K> &p, const Point_3<K> &q)
{
  return CGALi::z_equal(p, q, K());
}

CGAL_END_NAMESPACE

#endif  // CGAL_KERNEL_GLOBAL_FUNCTIONS_3_H
