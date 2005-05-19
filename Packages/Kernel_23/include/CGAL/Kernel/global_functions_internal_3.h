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
 
#ifndef CGAL_KERNEL_GLOBAL_FUNCTIONS_INTERNAL_3_H
#define CGAL_KERNEL_GLOBAL_FUNCTIONS_INTERNAL_3_H

// Generic functions calling the kernel functor.
// See comments in CGAL/Kernel/global_functions_internal_3.h.

CGAL_BEGIN_NAMESPACE

namespace CGALi {

template <typename K>
inline
Angle
angle(const typename CGAL_WRAP(K)::Point_3 &p,
      const typename CGAL_WRAP(K)::Point_3 &q,
      const typename CGAL_WRAP(K)::Point_3 &r, const K &k)
{
  return k.angle_3_object()(p, q, r);
}

template < class K >
inline
bool
are_ordered_along_line(const typename CGAL_WRAP(K)::Point_3 &p,
                       const typename CGAL_WRAP(K)::Point_3 &q,
                       const typename CGAL_WRAP(K)::Point_3 &r, const K& k)
{
  return k.are_ordered_along_line_3_object()(p, q, r);
}

template < class K >
inline
bool
are_strictly_ordered_along_line(const typename CGAL_WRAP(K)::Point_3 &p,
                                const typename CGAL_WRAP(K)::Point_3 &q,
                                const typename CGAL_WRAP(K)::Point_3 &r,
                                const K& k)
{
  return k.are_strictly_ordered_along_line_3_object()(p, q, r);
}

template <typename K>
inline
typename K::Plane_3
bisector(const typename CGAL_WRAP(K)::Point_3 &p,
         const typename CGAL_WRAP(K)::Point_3 &q, const K &k)
{
  return k.construct_bisector_3_object()(p, q);
}

template <typename K>
inline
typename K::Plane_3
bisector(const typename CGAL_WRAP(K)::Plane_3 &h1,
         const typename CGAL_WRAP(K)::Plane_3 &h2, const K &k)
{
  return k.construct_bisector_3_object()(h1, h2);
}

template < class K >
inline
typename K::Point_3
centroid(const typename CGAL_WRAP(K)::Point_3 &p,
         const typename CGAL_WRAP(K)::Point_3 &q,
         const typename CGAL_WRAP(K)::Point_3 &r,
         const typename CGAL_WRAP(K)::Point_3 &s, const K &k)
{
  return k.construct_centroid_3_object()(p, q, r, s);
}

template < class K >
inline
typename K::Point_3
centroid(const typename CGAL_WRAP(K)::Point_3 &p,
         const typename CGAL_WRAP(K)::Point_3 &q,
         const typename CGAL_WRAP(K)::Point_3 &r, const K &k)
{
  return k.construct_centroid_3_object()(p, q, r);
}

template < class K >
inline
typename K::Point_3
centroid(const typename CGAL_WRAP(K)::Tetrahedron_3 &t, const K &k)
{
  return k.construct_centroid_3_object()(t);
}

template < class K >
inline
typename K::Point_3
centroid(const typename CGAL_WRAP(K)::Triangle_3 &t, const K &k)
{
  return k.construct_centroid_3_object()(t);
}

template < class K >
inline
typename K::Point_3
circumcenter(const typename CGAL_WRAP(K)::Point_3 &p,
             const typename CGAL_WRAP(K)::Point_3 &q,
             const typename CGAL_WRAP(K)::Point_3 &r,
             const typename CGAL_WRAP(K)::Point_3 &s, const K &k)
{
  return k.construct_circumcenter_3_object()(p, q, r, s);
}

template < class K >
inline
typename K::Point_3
circumcenter(const typename CGAL_WRAP(K)::Tetrahedron_3 &t, const K& k)
{
  return k.construct_circumcenter_3_object()(t);
}

template < class K >
inline
typename K::Point_3
circumcenter(const typename CGAL_WRAP(K)::Point_3 &p,
             const typename CGAL_WRAP(K)::Point_3 &q,
             const typename CGAL_WRAP(K)::Point_3 &r, const K &k)
{
  return k.construct_circumcenter_3_object()(p, q, r);
}

template < class K >
inline
typename K::Point_3
circumcenter(const typename CGAL_WRAP(K)::Triangle_3 &t, const K& k)
{
  return k.construct_circumcenter_3_object()(t);
}

template < class K >
inline
bool
collinear(const typename CGAL_WRAP(K)::Point_3 &p,
          const typename CGAL_WRAP(K)::Point_3 &q,
          const typename CGAL_WRAP(K)::Point_3 &r,
          const K& k)
{
  return k.collinear_3_object()(p, q, r);
}

template < class K >
inline
bool
collinear_are_ordered_along_line(
          const typename CGAL_WRAP(K)::Point_3 &p,
          const typename CGAL_WRAP(K)::Point_3 &q,
          const typename CGAL_WRAP(K)::Point_3 &r,
          const K& k)
{
  return k.collinear_are_ordered_along_line_3_object()(p, q, r);
}

template < class K >
inline
bool
collinear_are_strictly_ordered_along_line(
          const typename CGAL_WRAP(K)::Point_3 &p,
          const typename CGAL_WRAP(K)::Point_3 &q,
          const typename CGAL_WRAP(K)::Point_3 &r,
          const K& k)
{
  return k.collinear_are_strictly_ordered_along_line_3_object()(p, q, r);
}

template < class K >
inline
Comparison_result
compare_distance_to_point(const typename CGAL_WRAP(K)::Point_3 &p,
                          const typename CGAL_WRAP(K)::Point_3 &q,
                          const typename CGAL_WRAP(K)::Point_3 &r,
			  const K& k)
{
  return k.compare_distance_3_object()(p, q, r);
}

template < class K >
inline
Comparison_result
compare_lexicographically_xyz(const typename CGAL_WRAP(K)::Point_3 &p,
                              const typename CGAL_WRAP(K)::Point_3 &q,
			      const K& k)
{
  return k.compare_xyz_3_object()(p, q);
}

template < class K >
inline
Comparison_result
compare_signed_distance_to_plane(const typename CGAL_WRAP(K)::Plane_3 &h,
				 const typename CGAL_WRAP(K)::Point_3 &p,
				 const typename CGAL_WRAP(K)::Point_3 &q,
				 const K &k)
{ 
  if (k.less_signed_distance_to_plane_3_object()(h, p, q)) return SMALLER;
  if (k.less_signed_distance_to_plane_3_object()(h, q, p)) return LARGER;
  return EQUAL;
}

template < class K >
inline
Comparison_result
compare_signed_distance_to_plane(const typename CGAL_WRAP(K)::Point_3 &hp,
				 const typename CGAL_WRAP(K)::Point_3 &hq,
				 const typename CGAL_WRAP(K)::Point_3 &hr,
				 const typename CGAL_WRAP(K)::Point_3 &p,
				 const typename CGAL_WRAP(K)::Point_3 &q,
				 const K &k)
{ 
  if (k.less_signed_distance_to_plane_3_object()(hp, hq, hr, p, q))
    return SMALLER;
  if (k.less_signed_distance_to_plane_3_object()(hp, hq, hr, q, p))
    return LARGER;
  return EQUAL;
}

template < class K >
inline
Comparison_result
compare_x(const typename CGAL_WRAP(K)::Point_3 &p,
	  const typename CGAL_WRAP(K)::Point_3 &q,
	  const K &k)
{ 
  return k.compare_x_3_object()(p, q);
}

template < class K >
inline
Comparison_result
compare_y(const typename CGAL_WRAP(K)::Point_3 &p,
	  const typename CGAL_WRAP(K)::Point_3 &q,
	  const K &k)
{ 
  return k.compare_y_3_object()(p, q);
}

template < class K >
inline
Comparison_result
compare_z(const typename CGAL_WRAP(K)::Point_3 &p,
	  const typename CGAL_WRAP(K)::Point_3 &q,
	  const K &k)
{ 
  return k.compare_z_3_object()(p, q);
}

template < class K >
inline
Comparison_result
compare_xyz(const typename CGAL_WRAP(K)::Point_3 &p,
	    const typename CGAL_WRAP(K)::Point_3 &q,
	    const K &k)
{ 
  return k.compare_xyz_3_object()(p, q);
}

template < class K >
inline
bool
coplanar(const typename CGAL_WRAP(K)::Point_3 &p,
         const typename CGAL_WRAP(K)::Point_3 &q,
         const typename CGAL_WRAP(K)::Point_3 &r,
         const typename CGAL_WRAP(K)::Point_3 &s, const K& k)
{
  return k.coplanar_3_object()(p, q, r, s);
}

template < class K >
inline
Orientation
coplanar_orientation(const typename CGAL_WRAP(K)::Point_3 &p,
                     const typename CGAL_WRAP(K)::Point_3 &q,
                     const typename CGAL_WRAP(K)::Point_3 &r,
                     const typename CGAL_WRAP(K)::Point_3 &s, const K& k)
{
  return k.coplanar_orientation_3_object()(p, q, r, s);
}

template < class K >
inline
Orientation
coplanar_orientation(const typename CGAL_WRAP(K)::Point_3 &p,
                     const typename CGAL_WRAP(K)::Point_3 &q,
                     const typename CGAL_WRAP(K)::Point_3 &r, const K& k)
{
  return k.coplanar_orientation_3_object()(p, q, r);
}

template < class K >
inline
Bounded_side
coplanar_side_of_bounded_circle(const typename CGAL_WRAP(K)::Point_3 &p,
                                const typename CGAL_WRAP(K)::Point_3 &q,
                                const typename CGAL_WRAP(K)::Point_3 &r,
                                const typename CGAL_WRAP(K)::Point_3 &t,
                                const K& k)
{
  return k.coplanar_side_of_bounded_circle_3_object()(p, q, r, t);
}

template < class K >
inline
typename K::Vector_3
cross_product(const typename CGAL_WRAP(K)::Vector_3 &v,
              const typename CGAL_WRAP(K)::Vector_3 &w, const K& k)
{
  return k.construct_cross_product_vector_3_object()(v, w);
}

template < class K >
inline
bool
has_smaller_distance_to_point(const typename CGAL_WRAP(K)::Point_3 &p,
                              const typename CGAL_WRAP(K)::Point_3 &q,
                              const typename CGAL_WRAP(K)::Point_3 &r,
			      const K &k)
{
  return k.less_distance_to_point_3_object()(p, q, r);
}

template < class K >
inline
bool
has_larger_distance_to_point(const typename CGAL_WRAP(K)::Point_3 &p,
			     const typename CGAL_WRAP(K)::Point_3 &q,
			     const typename CGAL_WRAP(K)::Point_3 &r,
			     const K &k)
{
  return k.compare_distance_3_object()(p, q, r) == LARGER;
}

template < class K >
inline
bool
has_larger_signed_distance_to_plane(const typename CGAL_WRAP(K)::Plane_3 &h,
				    const typename CGAL_WRAP(K)::Point_3 &p,
				    const typename CGAL_WRAP(K)::Point_3 &q,
				    const K &k)
{ 
  return k.less_signed_distance_to_plane_3_object()(h, q, p);
}

template < class K >
inline
bool
has_larger_signed_distance_to_plane(const typename CGAL_WRAP(K)::Point_3 &hp,
				    const typename CGAL_WRAP(K)::Point_3 &hq,
				    const typename CGAL_WRAP(K)::Point_3 &hr,
				    const typename CGAL_WRAP(K)::Point_3 &p,
				    const typename CGAL_WRAP(K)::Point_3 &q,
				    const K &k)
{ 
  return k.less_signed_distance_to_plane_3_object()(hp, hq, hr, q, p);
}

template < class K >
inline
bool
has_smaller_signed_distance_to_plane(const typename CGAL_WRAP(K)::Plane_3 &h,
                                     const typename CGAL_WRAP(K)::Point_3 &p,
                                     const typename CGAL_WRAP(K)::Point_3 &q,
				     const K &k)
{ 
  return k.less_signed_distance_to_plane_3_object()(h, p, q);
}

template < class K >
inline
bool
has_smaller_signed_distance_to_plane(const typename CGAL_WRAP(K)::Point_3 &hp,
                                     const typename CGAL_WRAP(K)::Point_3 &hq,
                                     const typename CGAL_WRAP(K)::Point_3 &hr,
                                     const typename CGAL_WRAP(K)::Point_3 &p,
                                     const typename CGAL_WRAP(K)::Point_3 &q,
				     const K &k)
{ 
  return k.less_signed_distance_to_plane_3_object()(hp, hq, hr, p, q);
}

template < class K >
inline
bool
less_x(const typename CGAL_WRAP(K)::Point_3 &p,
       const typename CGAL_WRAP(K)::Point_3 &q,
       const K &k)
{ 
  return k.less_x_3_object()(p, q);
}

template < class K >
inline
bool
less_y(const typename CGAL_WRAP(K)::Point_3 &p,
       const typename CGAL_WRAP(K)::Point_3 &q,
       const K &k)
{ 
  return k.less_y_3_object()(p, q);
}

template < class K >
inline
bool
less_z(const typename CGAL_WRAP(K)::Point_3 &p,
       const typename CGAL_WRAP(K)::Point_3 &q,
       const K &k)
{ 
  return k.less_z_3_object()(p, q);
}

template < class K >
inline
bool
lexicographically_xyz_smaller(const typename CGAL_WRAP(K)::Point_3 &p,
                              const typename CGAL_WRAP(K)::Point_3 &q,
                              const K &k)
{
  return k.less_xyz_3_object()(p, q);
}

template < class K >
inline
typename K::Point_3
midpoint(const typename CGAL_WRAP(K)::Point_3 &p,
         const typename CGAL_WRAP(K)::Point_3 &q, const K &k)
{
  return k.construct_midpoint_3_object()(p, q);
}

template < class K >
inline
typename K::Point_3
max_vertex(const typename CGAL_WRAP(K)::Iso_cuboid_3 &ic, const K &k)
{
  return k.construct_max_vertex_3_object()(ic);
}

template < class K >
inline
typename K::Point_3
min_vertex(const typename CGAL_WRAP(K)::Iso_cuboid_3 &ic, const K &k)
{
  return k.construct_min_vertex_3_object()(ic);
}

template <class K >
inline
Orientation
orientation(const typename CGAL_WRAP(K)::Point_3 &p,
	    const typename CGAL_WRAP(K)::Point_3 &q,
	    const typename CGAL_WRAP(K)::Point_3 &r,
	    const typename CGAL_WRAP(K)::Point_3 &s, const K &k)
{
  return k.orientation_3_object()(p, q, r, s);
}

template <class K >
inline
Orientation
orientation(const typename CGAL_WRAP(K)::Vector_3 &u,
	    const typename CGAL_WRAP(K)::Vector_3 &v,
	    const typename CGAL_WRAP(K)::Vector_3 &w, const K &k)
{
  return k.orientation_3_object()(u, v, w);
}

template < class K >
inline
typename K::Vector_3
orthogonal_vector(const typename CGAL_WRAP(K)::Point_3 &p,
		  const typename CGAL_WRAP(K)::Point_3 &q,
		  const typename CGAL_WRAP(K)::Point_3 &r, const K &k)
{
  return k.construct_orthogonal_vector_3_object()(p, q, r);
}

template < class K >
inline
typename K::Vector_3
orthogonal_vector(const typename CGAL_WRAP(K)::Plane_3 &p, const K &k)
{
  return k.construct_orthogonal_vector_3_object()(p);
}

template <typename K>
inline
bool
parallel(const typename CGAL_WRAP(K)::Line_3 &l1,
         const typename CGAL_WRAP(K)::Line_3 &l2, const K &k)
{
  return k.are_parallel_3_object()(l1, l2);
}

template <typename K>
inline
bool
parallel(const typename CGAL_WRAP(K)::Plane_3 &h1,
         const typename CGAL_WRAP(K)::Plane_3 &h2, const K &k)
{
  return k.are_parallel_3_object()(h1, h2);
}

template <typename K>
inline
bool
parallel(const typename CGAL_WRAP(K)::Ray_3 &r1,
         const typename CGAL_WRAP(K)::Ray_3 &r2, const K &k)
{
  return k.are_parallel_3_object()(r1, r2);
}

template <typename K>
inline
bool
parallel(const typename CGAL_WRAP(K)::Segment_3 &s1,
         const typename CGAL_WRAP(K)::Segment_3 &s2, const K &k)
{
  return k.are_parallel_3_object()(s1, s2);
}

template <class K >
inline
Bounded_side
side_of_bounded_sphere(const typename CGAL_WRAP(K)::Point_3 &p,
                       const typename CGAL_WRAP(K)::Point_3 &q,
                       const typename CGAL_WRAP(K)::Point_3 &test, const K &k)
{
  return k.side_of_bounded_sphere_3_object()(p, q, test);
}

template <class K >
inline
Bounded_side
side_of_bounded_sphere(const typename CGAL_WRAP(K)::Point_3 &p,
                       const typename CGAL_WRAP(K)::Point_3 &q,
                       const typename CGAL_WRAP(K)::Point_3 &r,
                       const typename CGAL_WRAP(K)::Point_3 &test, const K &k)
{
  return k.side_of_bounded_sphere_3_object()(p, q, r, test);
}

template <class K >
inline
Bounded_side
side_of_bounded_sphere(const typename CGAL_WRAP(K)::Point_3 &p,
                       const typename CGAL_WRAP(K)::Point_3 &q,
                       const typename CGAL_WRAP(K)::Point_3 &r,
                       const typename CGAL_WRAP(K)::Point_3 &s,
                       const typename CGAL_WRAP(K)::Point_3 &test, const K &k)
{
  return k.side_of_bounded_sphere_3_object()(p, q, r, s, test);
}

template <class K >
inline
Oriented_side
side_of_oriented_sphere(const typename CGAL_WRAP(K)::Point_3 &p,
                        const typename CGAL_WRAP(K)::Point_3 &q,
                        const typename CGAL_WRAP(K)::Point_3 &r,
                        const typename CGAL_WRAP(K)::Point_3 &s,
                        const typename CGAL_WRAP(K)::Point_3 &test, const K &k)
{
  return k.side_of_oriented_sphere_3_object()(p, q, r, s, test);
}

template <typename K>
inline
typename K::FT
squared_area(const typename CGAL_WRAP(K)::Point_3 &p,
	     const typename CGAL_WRAP(K)::Point_3 &q,
	     const typename CGAL_WRAP(K)::Point_3 &r, const K &k)
{
  return k.compute_squared_area_3_object()(p, q, r);
}

template < class K >
inline
typename K::FT
squared_radius(const typename CGAL_WRAP(K)::Point_3 &p,
	       const typename CGAL_WRAP(K)::Point_3 &q,
	       const typename CGAL_WRAP(K)::Point_3 &r,
	       const typename CGAL_WRAP(K)::Point_3 &s, const K &k)
{
  return k.compute_squared_radius_3_object()(p, q, r, s);
}

template < class K >
inline
typename K::FT
squared_radius(const typename CGAL_WRAP(K)::Point_3 &p,
	       const typename CGAL_WRAP(K)::Point_3 &q,
	       const typename CGAL_WRAP(K)::Point_3 &r, const K &k)
{
  return k.compute_squared_radius_3_object()(p, q, r);
}

template < class K >
inline
typename K::FT
squared_radius(const typename CGAL_WRAP(K)::Point_3 &p,
	       const typename CGAL_WRAP(K)::Point_3 &q, const K &k)
{
  return k.compute_squared_radius_3_object()(p, q);
}

template < class K >
inline
typename K::FT
volume(const typename CGAL_WRAP(K)::Point_3 &p,
       const typename CGAL_WRAP(K)::Point_3 &q,
       const typename CGAL_WRAP(K)::Point_3 &r,
       const typename CGAL_WRAP(K)::Point_3 &s, const K &k)
{
  return k.compute_volume_3_object()(p, q, r, s);
}

template < class K >
inline
bool
x_equal(const typename CGAL_WRAP(K)::Point_3 &p,
        const typename CGAL_WRAP(K)::Point_3 &q, const K &k)
{
  return k.equal_x_3_object()(p, q);
}

template < class K >
inline
bool
y_equal(const typename CGAL_WRAP(K)::Point_3 &p,
        const typename CGAL_WRAP(K)::Point_3 &q, const K &k)
{
  return k.equal_y_3_object()(p, q);
}

template < class K >
inline
bool
z_equal(const typename CGAL_WRAP(K)::Point_3 &p,
        const typename CGAL_WRAP(K)::Point_3 &q, const K &k)
{
  return k.equal_z_3_object()(p, q);
}


// The following functions only call some of the previous ones.

template <typename K>
inline
bool
are_negative_oriented(const typename CGAL_WRAP(K)::Point_3 &p,
                      const typename CGAL_WRAP(K)::Point_3 &q,
                      const typename CGAL_WRAP(K)::Point_3 &r,
                      const typename CGAL_WRAP(K)::Point_3 &s, const K &k)
{
  return CGALi::orientation(p, q, r, s, k) == NEGATIVE;
}

template <typename K>
inline
bool
are_positive_oriented(const typename CGAL_WRAP(K)::Point_3 &p,
                      const typename CGAL_WRAP(K)::Point_3 &q,
                      const typename CGAL_WRAP(K)::Point_3 &r,
                      const typename CGAL_WRAP(K)::Point_3 &s, const K &k)
{
  return CGALi::orientation(p, q, r, s, k) == POSITIVE;
}

template < class K >
inline
bool
lexicographically_xyz_smaller_or_equal(const typename CGAL_WRAP(K)::Point_3 &p,
                                       const typename CGAL_WRAP(K)::Point_3 &q,
                                       const K&k)
{
  return CGALi::compare_lexicographically_xyz(p, q, k) != LARGER;
}

} // namespace CGALi

CGAL_END_NAMESPACE

#endif  // CGAL_KERNEL_GLOBAL_FUNCTIONS_INTERNAL_3_H
