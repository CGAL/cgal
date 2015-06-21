// Copyright (c) 2003-2004  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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

namespace CGAL {

template <typename K>
inline
Angle
angle(const Vector_3<K> &u, const Vector_3<K> &v)
{
  return internal::angle(u, v, K());
}

template <typename K>
inline
Angle
angle(const Point_3<K> &p, const Point_3<K> &q, const Point_3<K> &r)
{
  return internal::angle(p, q, r, K());
}

template <typename K>
inline
Angle
angle(const Point_3<K> &p, const Point_3<K> &q,
      const Point_3<K> &r, const Point_3<K> &s)
{
  return internal::angle(p, q, r, s, K());
}

template < typename K >
inline
typename K::Boolean
are_negative_oriented(const Point_3<K>& p, const Point_3<K>& q,
                      const Point_3<K>& r, const Point_3<K>& s)
{
  return internal::are_negative_oriented(p, q, r, s, K());
}

template < class K >
inline
typename K::Boolean
are_ordered_along_line(const Point_3<K> &p,
                       const Point_3<K> &q,
                       const Point_3<K> &r)
{
  return internal::are_ordered_along_line(p, q, r, K());
}

template < typename K >
inline
typename K::Boolean
are_positive_oriented(const Point_3<K>& p, const Point_3<K>& q,
                      const Point_3<K>& r, const Point_3<K>& s)
{
  return internal::are_positive_oriented(p, q, r, s, K());
}

template < class K >
inline
typename K::Boolean
are_strictly_ordered_along_line(const Point_3<K> &p,
                                const Point_3<K> &q,
                                const Point_3<K> &r)
{
  return internal::are_strictly_ordered_along_line(p, q, r, K());
}

template < class K >
inline
typename K::Point_3
barycenter(const Point_3<K> &p1, const typename K::FT& w1,
           const Point_3<K> &p2)
{
  return internal::barycenter(p1, w1, p2, K());
}

template < class K >
inline
typename K::Point_3
barycenter(const Point_3<K> &p1, const typename K::FT& w1,
           const Point_3<K> &p2, const typename K::FT& w2)
{
  return internal::barycenter(p1, w1, p2, w2, K());
}

template < class K >
inline
typename K::Point_3
barycenter(const Point_3<K> &p1, const typename K::FT& w1,
           const Point_3<K> &p2, const typename K::FT& w2,
           const Point_3<K> &p3)
{
  return internal::barycenter(p1, w1, p2, w2, p3, K());
}

template < class K >
inline
typename K::Point_3
barycenter(const Point_3<K> &p1, const typename K::FT& w1,
           const Point_3<K> &p2, const typename K::FT& w2,
           const Point_3<K> &p3, const typename K::FT& w3)
{
  return internal::barycenter(p1, w1, p2, w2, p3, w3, K());
}

template < class K >
inline
typename K::Point_3
barycenter(const Point_3<K> &p1, const typename K::FT& w1,
           const Point_3<K> &p2, const typename K::FT& w2,
           const Point_3<K> &p3, const typename K::FT& w3,
           const Point_3<K> &p4)
{
  return internal::barycenter(p1, w1, p2, w2, p3, w3, p4, K());
}

template < class K >
inline
typename K::Point_3
barycenter(const Point_3<K> &p1, const typename K::FT& w1,
           const Point_3<K> &p2, const typename K::FT& w2,
           const Point_3<K> &p3, const typename K::FT& w3,
           const Point_3<K> &p4, const typename K::FT& w4)
{
  return internal::barycenter(p1, w1, p2, w2, p3, w3, p4, w4, K());
}

template <typename K>
inline
typename K::Plane_3
bisector(const Point_3<K> &p, const Point_3<K> &q)
{
  return internal::bisector(p, q, K());
}

template <typename K>
inline
typename K::Plane_3
bisector(const Plane_3<K> &h1, const Plane_3<K> &h2)
{
  return internal::bisector(h1, h2, K());
}

template < class K >
inline
Point_3<K>
centroid(const Point_3<K> &p, const Point_3<K> &q,
         const Point_3<K> &r, const Point_3<K> &s)
{
  return internal::centroid(p, q, r, s, K());
}

template < class K >
inline
Point_3<K>
centroid(const Point_3<K> &p, const Point_3<K> &q, const Point_3<K> &r)
{
  return internal::centroid(p, q, r, K());
}

template < class K >
inline
Point_3<K>
centroid(const Tetrahedron_3<K> &t)
{
  return internal::centroid(t, K());
}

template < class K >
inline
Point_3<K>
centroid(const Triangle_3<K> &t)
{
  return internal::centroid(t, K());
}

template < class K >
inline
typename K::Point_3
circumcenter(const Point_3<K> &p,
             const Point_3<K> &q)
{
  return internal::circumcenter(p, q, K());
}

template < class K >
inline
typename K::Point_3
circumcenter(const Point_3<K> &p,
             const Point_3<K> &q,
             const Point_3<K> &r)
{
  return internal::circumcenter(p, q, r, K());
}

template < class K >
inline
typename K::Point_3
circumcenter(const Point_3<K> &p, const Point_3<K> &q,
             const Point_3<K> &r, const Point_3<K> &s)
{
  return internal::circumcenter(p, q, r, s, K());
}

template < class K >
inline
typename K::Point_3
circumcenter(const Tetrahedron_3<K> &t)
{
  return internal::circumcenter(t, K());
}

template < class K >
inline
typename K::Point_3
circumcenter(const Triangle_3<K> &t)
{
  return internal::circumcenter(t, K());
}

template < class K >
inline
typename K::Boolean
collinear(const Point_3<K> &p, const Point_3<K> &q, const Point_3<K> &r)
{
  return internal::collinear(p, q, r, K());
}

template < class K >
inline
typename K::Boolean
collinear_are_ordered_along_line(const Point_3<K> &p,
                                 const Point_3<K> &q,
                                 const Point_3<K> &r)
{
  return internal::collinear_are_ordered_along_line(p, q, r, K());
}

template < class K >
inline
typename K::Boolean
collinear_are_strictly_ordered_along_line(const Point_3<K> &p,
                                          const Point_3<K> &q,
                                          const Point_3<K> &r)
{
  return internal::collinear_are_strictly_ordered_along_line(p, q, r, K());
}

template < class K >
inline
typename K::Comparison_result
compare_dihedral_angle(const Point_3<K>& a1, const Point_3<K>& b1, 
                       const Point_3<K>& c1, const Point_3<K>& d1, 
                       const Point_3<K>& a2, const Point_3<K>& b2, 
                       const Point_3<K>& c2, const Point_3<K>& d2)
{
  return internal::compare_dihedral_angle(a1, b1, c1, d1, a2, b2, c2, d2, K());
}

template < class K >
inline
typename K::Comparison_result
compare_dihedral_angle(const Point_3<K>& a1, const Point_3<K>& b1, 
                       const Point_3<K>& c1, const Point_3<K>& d1, 
                       const typename K::FT& cosine)
{
  return internal::compare_dihedral_angle(a1, b1, c1, d1, cosine, K());
}

template < class K >
inline
typename K::Comparison_result
compare_dihedral_angle(const Vector_3<K>& ab1, 
                       const Vector_3<K>& ac1,
                       const Vector_3<K>& ad1,
                       const Vector_3<K>& ab2,
                       const Vector_3<K>& ac2,
                       const Vector_3<K>& ad2)
{
  return internal::compare_dihedral_angle(ab1, ac1, ad1, ab2, ac2, ad2, K());
}

template < class K >
inline
typename K::Comparison_result
compare_dihedral_angle(const Vector_3<K>& ab1, 
                       const Vector_3<K>& ac1,
                       const Vector_3<K>& ad1,
                       const typename K::FT& cosine)
{
  return internal::compare_dihedral_angle(ab1, ac1, ad1, cosine, K());
}

template < class K >
inline
typename K::Comparison_result
compare_distance_to_point(const Point_3<K> &p,
                          const Point_3<K> &q,
                          const Point_3<K> &r)
{
  return internal::compare_distance_to_point(p, q, r, K());
}

template < class K >
inline
typename K::Comparison_result
compare_squared_distance(const Point_3<K> &p,
                         const Point_3<K> &q,
                         const typename K::FT &d2)
{
  return internal::compare_squared_distance(p, q, d2, K());
}

template < class K >
inline
typename K::Comparison_result
compare_squared_radius(const Point_3<K> &p,
		       const typename K::FT &sr)
{
  return internal::compare_squared_radius(p, sr, K());
}

template < class K >
inline
typename K::Comparison_result
compare_squared_radius(const Point_3<K> &p,
		       const Point_3<K> &q,
		       const typename K::FT &sr)
{
  return internal::compare_squared_radius(p, q, sr, K());
}

template < class K >
inline
typename K::Comparison_result
compare_squared_radius(const Point_3<K> &p,
		       const Point_3<K> &q,
		       const Point_3<K> &r,
		       const typename K::FT &sr)
{
  return internal::compare_squared_radius(p, q, r, sr, K());
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
  return internal::compare_squared_radius(p, q, r, s, sr, K());
}

template < class K >
inline
typename K::Comparison_result
compare_lexicographically_xyz(const Point_3<K> &p,
                              const Point_3<K> &q)
{
  return internal::compare_lexicographically_xyz(p, q, K());
}

template < class K >
inline
typename K::Comparison_result
compare_lexicographically(const Point_3<K> &p,
                          const Point_3<K> &q)
{
    return internal::compare_lexicographically_xyz(p, q, K());
}

template < class K >
inline
typename K::Comparison_result
compare_signed_distance_to_plane(const Plane_3<K> &h,
				 const Point_3<K> &p,
				 const Point_3<K> &q)
{ 
  return internal::compare_signed_distance_to_plane(h, p, q, K());
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
  return internal::compare_signed_distance_to_plane(hp, hq, hr, p, q, K());
}

template < class K >
inline
typename K::Comparison_result
compare_x(const Point_3<K> &p, const Point_3<K> &q)
{ 
  return internal::compare_x(p, q, K());
}

template < class K >
inline
typename K::Comparison_result
compare_y(const Point_3<K> &p, const Point_3<K> &q)
{ 
  return internal::compare_y(p, q, K());
}

template < class K >
inline
typename K::Comparison_result
compare_z(const Point_3<K> &p, const Point_3<K> &q)
{ 
  return internal::compare_z(p, q, K());
}

template < class K >
inline
typename K::Comparison_result
compare_xyz(const Point_3<K> &p, const Point_3<K> &q)
{ 
  return internal::compare_xyz(p, q, K());
}

template < class K >
inline
typename K::Boolean
coplanar(const Point_3<K> &p, const Point_3<K> &q,
         const Point_3<K> &r, const Point_3<K> &s)
{
  return internal::coplanar(p, q, r, s, K());
}


template < class K >
inline
typename K::Orientation
coplanar_orientation(const Point_3<K> &p,
                     const Point_3<K> &q,
                     const Point_3<K> &r,
                     const Point_3<K> &s)
{
  return internal::coplanar_orientation(p, q, r, s, K());
}

template < class K >
inline
typename K::Orientation
coplanar_orientation(const Point_3<K> &p,
                     const Point_3<K> &q,
                     const Point_3<K> &r)
{
  return internal::coplanar_orientation(p, q, r, K());
}

template < class K >
inline
typename K::Bounded_side
coplanar_side_of_bounded_circle(const Point_3<K> &p,
                                const Point_3<K> &q,
                                const Point_3<K> &r,
                                const Point_3<K> &t)
{
  return internal::coplanar_side_of_bounded_circle(p, q, r, t, K());
}

template < class K >
inline
typename K::Vector_3
cross_product(const Vector_3<K> &v, const Vector_3<K> &w)
{
  return internal::cross_product(v, w, K());
}

template < class K >
inline
typename K::FT
determinant(const Vector_3<K> &v0, const Vector_3<K> &v1,
            const Vector_3<K> &v2)
{
  return internal::determinant(v0, v1, v2, K());
}

template < class K >
inline
typename K::Line_3
equidistant_line(const Point_3<K> &p, const Point_3<K> &q, const Point_3<K> &r)
{
  return internal::equidistant_line(p, q, r, K());
}

template < class K >
inline
typename K::Boolean
has_larger_distance_to_point(const Point_3<K> &p,
			     const Point_3<K> &q,
			     const Point_3<K> &r)
{
  return internal::has_larger_distance_to_point(p, q, r, K());
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
  return internal::has_larger_signed_distance_to_plane(hp, hq, hr, p, q, K());
}

template < class K >
inline
typename K::Boolean
has_larger_signed_distance_to_plane(const Plane_3<K> &h,
				    const Point_3<K> &p,
				    const Point_3<K> &q)
{ 
  return internal::has_larger_signed_distance_to_plane(h, p, q, K());
}

template < class K >
inline
typename K::Boolean
has_smaller_distance_to_point(const Point_3<K> &p,
                              const Point_3<K> &q,
                              const Point_3<K> &r)
{
  return internal::has_smaller_distance_to_point(p, q, r, K());
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
  return internal::has_smaller_signed_distance_to_plane(hp, hq, hr, p, q, K());
}

template < class K >
inline
typename K::Boolean
has_smaller_signed_distance_to_plane(const Plane_3<K> &h,
                                     const Point_3<K> &p,
                                     const Point_3<K> &q)
{ 
  return internal::has_smaller_signed_distance_to_plane(h, p, q, K());
}

template < class K >
inline
typename K::Boolean
less_x(const Point_3<K> &p, const Point_3<K> &q)
{ 
  return internal::less_x(p, q, K());
}

template < class K >
inline
typename K::Boolean
less_y(const Point_3<K> &p, const Point_3<K> &q)
{ 
  return internal::less_y(p, q, K());
}

template < class K >
inline
typename K::Boolean
less_z(const Point_3<K> &p, const Point_3<K> &q)
{ 
  return internal::less_z(p, q, K());
}

template < class K >
inline
typename K::Boolean
lexicographically_xyz_smaller(const Point_3<K> &p, const Point_3<K> &q)
{
  return internal::lexicographically_xyz_smaller(p, q, K());
}

template < class K >
inline
typename K::Boolean
lexicographically_xyz_smaller_or_equal(const Point_3<K> &p,
                                       const Point_3<K> &q)
{
  return internal::lexicographically_xyz_smaller_or_equal(p, q, K());
}

template < class K >
inline
typename K::Point_3
midpoint(const Point_3<K> &p, const Point_3<K> &q)
{
  return internal::midpoint(p, q, K());
}

template < class K >
inline
typename K::Point_3
max_vertex(const Iso_cuboid_3<K> &ic)
{
  return internal::max_vertex(ic, K());
}

template < class K >
inline
typename K::Point_3
min_vertex(const Iso_cuboid_3<K> &ic)
{
  return internal::min_vertex(ic, K());
}

template < class K >
inline
typename K::Vector_3
normal(const Point_3<K> &p, const Point_3<K> &q, const Point_3<K> &r)
{
  return internal::normal(p, q, r, K());
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
  return internal::orientation(p, q, r, s, K());
}

template <class K >
inline
typename K::Orientation
orientation(const Vector_3<K> &u, const Vector_3<K> &v, const Vector_3<K> &w)
{
  return internal::orientation(u, v, w, K());
}

template <class K >
inline
typename K::Vector_3
orthogonal_vector(const Point_3<K>& p,
		  const Point_3<K>& q,
		  const Point_3<K>& r)
{
  return internal::orthogonal_vector(p, q, r, K());
}

template <class K >
inline
typename K::Vector_3
orthogonal_vector(const Plane_3<K>& p)
{
  return internal::orthogonal_vector(p, K());
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
  return internal::side_of_bounded_sphere(p, q, test, K());
}

template <class K >
inline
typename K::Bounded_side
side_of_bounded_sphere(const Point_3<K> &p,
                       const Point_3<K> &q,
                       const Point_3<K> &r,
                       const Point_3<K> &test)
{
  return internal::side_of_bounded_sphere(p, q, r, test, K());
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
  return internal::side_of_bounded_sphere(p, q, r, s, test, K());
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
  return internal::side_of_oriented_sphere(p, q, r, s, test, K());
}

template <typename K>
inline
typename K::FT
squared_area(const Point_3<K> &p, const Point_3<K> &q, const Point_3<K> &r)
{
  return internal::squared_area(p, q, r, K());
}

template < class K >
inline
typename K::FT
squared_radius(const Point_3<K> &p, const Point_3<K> &q,
	       const Point_3<K> &r, const Point_3<K> &s)
{
  return internal::squared_radius(p, q, r, s, K());
}

template < class K >
inline
typename K::FT
squared_radius(const Point_3<K> &p, const Point_3<K> &q, const Point_3<K> &r)
{
  return internal::squared_radius(p, q, r, K());
}

template < class K >
inline
typename K::FT
squared_radius(const Point_3<K> &p, const Point_3<K> &q)
{
  return internal::squared_radius(p, q, K());
}

template < class K >
inline
typename K::FT
squared_radius(const Point_3<K> &p)
{
  return internal::squared_radius(p, K());
}

template < class K >
inline
typename K::Vector_3
unit_normal(const Point_3<K> &p, const Point_3<K> &q, const Point_3<K> &r)
{
  return internal::unit_normal(p, q, r, K());
}

template < class K >
inline
typename K::FT
volume(const Point_3<K> &p, const Point_3<K> &q,
       const Point_3<K> &r, const Point_3<K> &s)
{
  return internal::volume(p, q, r, s, K());
}

template < class K >
inline
typename K::Boolean
x_equal(const Point_3<K> &p, const Point_3<K> &q)
{
  return internal::x_equal(p, q, K());
}

template < class K >
inline
typename K::Boolean
y_equal(const Point_3<K> &p, const Point_3<K> &q)
{
  return internal::y_equal(p, q, K());
}

template < class K >
inline
typename K::Boolean
z_equal(const Point_3<K> &p, const Point_3<K> &q)
{
  return internal::z_equal(p, q, K());
}

} //namespace CGAL

#endif  // CGAL_KERNEL_GLOBAL_FUNCTIONS_3_H
