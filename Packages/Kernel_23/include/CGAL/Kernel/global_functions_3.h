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
// $Source$
// $Revision$ $Date$
// $Name$
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
bool
are_negative_oriented(const Point_3<K>& p, const Point_3<K>& q,
                      const Point_3<K>& r, const Point_3<K>& s)
{
  return CGALi::are_negative_oriented(p, q, r, s, K());
}

template < class K >
inline
bool
are_ordered_along_line(const Point_3<K> &p,
                       const Point_3<K> &q,
                       const Point_3<K> &r)
{
  return CGALi::are_ordered_along_line(p, q, r, K());
}

template < typename K >
inline
bool
are_positive_oriented(const Point_3<K>& p, const Point_3<K>& q,
                      const Point_3<K>& r, const Point_3<K>& s)
{
  return CGALi::are_positive_oriented(p, q, r, s, K());
}

template < class K >
inline
bool
are_strictly_ordered_along_line(const Point_3<K> &p,
                                const Point_3<K> &q,
                                const Point_3<K> &r)
{
  return CGALi::are_strictly_ordered_along_line(p, q, r, K());
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
circumcenter(const Point_3<K> &p,
             const Point_3<K> &q,
             const Point_3<K> &r)
{
  return CGALi::circumcenter(p, q, r, K());
}

template < class K >
inline
typename K::Point_3
circumcenter(const Triangle_3<K> &t)
{
  return CGALi::circumcenter(t, K());
}

template < class K >
CGAL_KERNEL_INLINE
bool
collinear(const Point_3<K> &p, const Point_3<K> &q, const Point_3<K> &r)
{
  return CGALi::collinear(p, q, r, K());
}

template < class K >
inline
Comparison_result
compare_distance_to_point(const Point_3<K> &p,
                          const Point_3<K> &q,
                          const Point_3<K> &r)
{
  return CGALi::compare_distance_to_point(p, q, r, K());
}

template < class K >
inline
Comparison_result
compare_signed_distance_to_plane(const Plane_3<K> &h,
				 const Point_3<K> &p,
				 const Point_3<K> &q)
{ 
  return CGALi::compare_signed_distance_to_plane(h, p, q, K());
}

template < class K >
inline
Comparison_result
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
bool
collinear_are_strictly_ordered_along_line(const Point_3<K> &p,
                                          const Point_3<K> &q,
                                          const Point_3<K> &r)
{
  return CGALi::collinear_are_strictly_ordered_along_line(p, q, r, K());
}

template < class K >
inline
bool
coplanar(const Point_3<K> &p, const Point_3<K> &q,
         const Point_3<K> &r, const Point_3<K> &s)
{
  return CGALi::coplanar(p, q, r, s, K());
}

template < class K >
inline
bool
has_larger_distance_to_point(const Point_3<K> &p,
			     const Point_3<K> &q,
			     const Point_3<K> &r)
{
  return CGALi::has_larger_distance_to_point(p, q, r, K());
}

template < class K >
inline
bool
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
bool
has_larger_signed_distance_to_plane(const Plane_3<K> &h,
				    const Point_3<K> &p,
				    const Point_3<K> &q)
{ 
  return CGALi::has_larger_signed_distance_to_plane(h, p, q, K());
}

template < class K >
inline
bool
has_smaller_distance_to_point(const Point_3<K> &p,
                              const Point_3<K> &q,
                              const Point_3<K> &r)
{
  return CGALi::has_smaller_distance_to_point(p, q, r, K());
}

template < class K >
inline
bool
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
bool
has_smaller_signed_distance_to_plane(const Plane_3<K> &h,
                                     const Point_3<K> &p,
                                     const Point_3<K> &q)
{ 
  return CGALi::has_smaller_signed_distance_to_plane(h, p, q, K());
}


template < class K >
inline
Orientation
coplanar_orientation(const Point_3<K> &p,
                     const Point_3<K> &q,
                     const Point_3<K> &r,
                     const Point_3<K> &s)
{
  return CGALi::coplanar_orientation(p, q, r, s, K());
}

template < class K >
inline
Orientation
coplanar_orientation(const Point_3<K> &p,
                     const Point_3<K> &q,
                     const Point_3<K> &r)
{
  return CGALi::coplanar_orientation(p, q, r, K());
}

template < class K >
inline
Bounded_side
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

// FIXME TODO : what to do with teh operators ?
template < class K >
inline
bool
operator<(const Point_3<K>& p, const Point_3<K>& q)
{ return K().less_xyz_3_object()(p, q); }

template < class K >
inline
bool
operator>(const Point_3<K>& p, const Point_3<K>& q)
{ return K().less_xyz_3_object()(q, p); }

template < class K >
inline
bool
operator<=(const Point_3<K>& p, const Point_3<K>& q)
{ return ! K().less_xyz_3_object()(q, p); }

template < class K >
inline
bool
operator>=(const Point_3<K>& p, const Point_3<K>& q)
{ return ! K().less_xyz_3_object()(p, q); }

template < class K >
inline
typename K::Vector_3
operator*(const typename CGAL_WRAP(K)::FT &c, const Vector_3<K> &w)
{
  return K().construct_scaled_vector_3_object()(w, c);
}

template < class K >
inline
typename K::Vector_3
operator*(const Vector_3<K> &w, const typename CGAL_WRAP(K)::FT &c)
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
Orientation
orientation(const Point_3<K> &p,
            const Point_3<K> &q,
            const Point_3<K> &r,
            const Point_3<K> &s)
{
  return CGALi::orientation(p, q, r, s, K());
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

template <class K >
inline
Bounded_side
side_of_bounded_sphere(const Point_3<K> &p,
                       const Point_3<K> &q,
                       const Point_3<K> &test)
{
  return CGALi::side_of_bounded_sphere(p, q, test, K());
}

template <class K >
inline
Bounded_side
side_of_bounded_sphere(const Point_3<K> &p,
                       const Point_3<K> &q,
                       const Point_3<K> &r,
                       const Point_3<K> &test)
{
  return CGALi::side_of_bounded_sphere(p, q, r, test, K());
}

template <class K >
inline
Bounded_side
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
Oriented_side
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
volume(const Point_3<K> &p, const Point_3<K> &q,
       const Point_3<K> &r, const Point_3<K> &s)
{
  return CGALi::volume(p, q, r, s, K());
}

CGAL_END_NAMESPACE

#endif  // CGAL_KERNEL_GLOBAL_FUNCTIONS_3_H
