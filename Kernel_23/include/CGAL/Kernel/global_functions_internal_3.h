// Copyright (c) 2003-2004  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// Copyright (c) 2010 GeometryFactory Sarl (France) 
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
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Sylvain Pion
 
#ifndef CGAL_KERNEL_GLOBAL_FUNCTIONS_INTERNAL_3_H
#define CGAL_KERNEL_GLOBAL_FUNCTIONS_INTERNAL_3_H

// Generic functions calling the kernel functor.
// See comments in CGAL/Kernel/global_functions_internal_3.h.

#include <CGAL/basic.h>
#include <CGAL/Dimension.h>
#include <boost/utility/enable_if.hpp>
#include <boost/mpl/equal_to.hpp>
#include <boost/mpl/integral_c.hpp>

namespace CGAL {

namespace internal {

template <typename K>
inline
typename K::Angle
angle(const typename K::Vector_3 &u,
      const typename K::Vector_3 &v,
      const K &k)
{
  return k.angle_3_object()(u, v);
}

template <typename K>
inline
typename K::Angle
angle(const typename K::Point_3 &p,
      const typename K::Point_3 &q,
      const typename K::Point_3 &r, const K &k)
{
  return k.angle_3_object()(p, q, r);
}

template <typename K>
inline
typename K::Angle
angle(const typename K::Point_3 &p,
      const typename K::Point_3 &q,
      const typename K::Point_3 &r,
      const typename K::Point_3 &s,
      const K &k)
{
  return k.angle_3_object()(p, q, r, s);
}

template <typename K>
inline
typename K::Angle
angle(const typename K::Point_3 &p,
      const typename K::Point_3 &q,
      const typename K::Point_3 &r,
      const typename K::Vector_3 &v,
      const K &k)
{
  return k.angle_3_object()(p, q, r, v);
}

template < class K >
inline
typename K::FT
approximate_dihedral_angle(const typename K::Point_3 &p,
                           const typename K::Point_3 &q,
                           const typename K::Point_3 &r,
                           const typename K::Point_3 &s, const K& k)
{
  return k.compute_approximate_dihedral_angle_3_object()(p, q, r, s);
}

template < class K >
inline
typename K::Boolean
are_ordered_along_line(const typename K::Point_3 &p,
                       const typename K::Point_3 &q,
                       const typename K::Point_3 &r, const K& k)
{
  return k.are_ordered_along_line_3_object()(p, q, r);
}

template < class K >
inline
typename K::Boolean
are_strictly_ordered_along_line(const typename K::Point_3 &p,
                                const typename K::Point_3 &q,
                                const typename K::Point_3 &r,
                                const K& k)
{
  return k.are_strictly_ordered_along_line_3_object()(p, q, r);
}

template < class K >
inline
typename K::Point_3
barycenter(const typename K::Point_3 &p1, const typename K::FT& w1,
           const typename K::Point_3 &p2, const K& k)
{
  return k.construct_barycenter_3_object()(p1, w1, p2);
}
  
template < class K >
inline
typename K::Point_3
barycenter(const typename K::Point_3 &p1, const typename K::FT& w1,
           const typename K::Point_3 &p2, const typename K::FT& w2, const K& k)
{
  return k.construct_barycenter_3_object()(p1, w1, p2, w2);
}
  
template < class K >
inline
typename K::Point_3
barycenter(const typename K::Point_3 &p1, const typename K::FT& w1,
           const typename K::Point_3 &p2, const typename K::FT& w2,
           const typename K::Point_3 &p3, const K& k)
{        
  return k.construct_barycenter_3_object()(p1, w1, p2, w2, p3);
}
  
template < class K >
inline
typename K::Point_3
barycenter(const typename K::Point_3 &p1, const typename K::FT& w1,
           const typename K::Point_3 &p2, const typename K::FT& w2,
           const typename K::Point_3 &p3, const typename K::FT& w3, const K& k)
{        
  return k.construct_barycenter_3_object()(p1, w1, p2, w2, p3, w3);
}        

template < class K >
inline
typename K::Point_3
barycenter(const typename K::Point_3 &p1, const typename K::FT& w1,
           const typename K::Point_3 &p2, const typename K::FT& w2,
           const typename K::Point_3 &p3, const typename K::FT& w3,
           const typename K::Point_3 &p4, const K& k)
{
  return k.construct_barycenter_3_object()(p1, w1, p2, w2, p3, w3, p4);
}

template < class K >
inline
typename K::Point_3
barycenter(const typename K::Point_3 &p1, const typename K::FT& w1,
           const typename K::Point_3 &p2, const typename K::FT& w2,
           const typename K::Point_3 &p3, const typename K::FT& w3,
           const typename K::Point_3 &p4, const typename K::FT& w4, const K& k)
{
  return k.construct_barycenter_3_object()(p1, w1, p2, w2, p3, w3, p4, w4);
}

template <typename K>
inline
typename K::Plane_3
bisector(const typename K::Point_3 &p,
         const typename K::Point_3 &q, const K &k)
{
  return k.construct_bisector_3_object()(p, q);
}

template <typename K>
inline
typename K::Plane_3
bisector(const typename K::Plane_3 &h1,
         const typename K::Plane_3 &h2, const K &k)
{
  return k.construct_bisector_3_object()(h1, h2);
}

template < class K >
inline
typename K::Point_3
centroid(const typename K::Point_3 &p,
         const typename K::Point_3 &q,
         const typename K::Point_3 &r,
         const typename K::Point_3 &s, const K &k)
{
  return k.construct_centroid_3_object()(p, q, r, s);
}

template < class K >
inline
typename K::Point_3
centroid(const typename K::Point_3 &p,
         const typename K::Point_3 &q,
         const typename K::Point_3 &r, const K &k)
{
  return k.construct_centroid_3_object()(p, q, r);
}

template < class K >
inline
typename K::Point_3
centroid(const typename K::Tetrahedron_3 &t, const K &k)
{
  return k.construct_centroid_3_object()(t);
}

template < class K >
inline
typename K::Point_3
centroid(const typename K::Triangle_3 &t, const K &k)
{
  return k.construct_centroid_3_object()(t);
}

template < class K >
inline
typename K::Point_3
circumcenter(const typename K::Point_3 &p,
             const typename K::Point_3 &q, const K &k)
{
  return k.construct_circumcenter_3_object()(p, q);
}

template < class K >
inline
typename K::Point_3
circumcenter(const typename K::Point_3 &p,
             const typename K::Point_3 &q,
             const typename K::Point_3 &r, const K &k)
{
  return k.construct_circumcenter_3_object()(p, q, r);
}

template < class K >
inline
typename K::Point_3
circumcenter(const typename K::Point_3 &p,
             const typename K::Point_3 &q,
             const typename K::Point_3 &r,
             const typename K::Point_3 &s, const K &k)
{
  return k.construct_circumcenter_3_object()(p, q, r, s);
}

template < class K >
inline
typename K::Point_3
circumcenter(const typename K::Tetrahedron_3 &t, const K& k)
{
  return k.construct_circumcenter_3_object()(t);
}

template < class K >
inline
typename K::Point_3
circumcenter(const typename K::Triangle_3 &t, const K& k)
{
  return k.construct_circumcenter_3_object()(t);
}

template < class K >
inline
typename K::Boolean
collinear(const typename K::Point_3 &p,
          const typename K::Point_3 &q,
          const typename K::Point_3 &r,
          const K& k)
{
  return k.collinear_3_object()(p, q, r);
}

template < class K >
inline
typename K::Boolean
collinear_are_ordered_along_line(
          const typename K::Point_3 &p,
          const typename K::Point_3 &q,
          const typename K::Point_3 &r,
          const K& k)
{
  return k.collinear_are_ordered_along_line_3_object()(p, q, r);
}

template < class K >
inline
typename K::Boolean
collinear_are_strictly_ordered_along_line(
          const typename K::Point_3 &p,
          const typename K::Point_3 &q,
          const typename K::Point_3 &r,
          const K& k)
{
  return k.collinear_are_strictly_ordered_along_line_3_object()(p, q, r);
}


template < class K >
inline
typename K::Comparison_result
compare_dihedral_angle(const typename K::Point_3& a1,
                       const typename K::Point_3& b1, 
                       const typename K::Point_3& c1,
                       const typename K::Point_3& d1, 
                       const typename K::Point_3& a2, 
                       const typename K::Point_3& b2, 
                       const typename K::Point_3& c2,
                       const typename K::Point_3& d2,
                       const K& k)
{
  return k.compare_dihedral_angle_3_object()(a1, b1, c1, d1, a2, b2, c2, d2);
}

template < class K >
inline
typename K::Comparison_result
compare_dihedral_angle(const typename K::Point_3& a1,
                       const typename K::Point_3& b1, 
                       const typename K::Point_3& c1,
                       const typename K::Point_3& d1, 
                       const typename K::FT& cosine,
                       const K& k)
{
  return k.compare_dihedral_angle_3_object()(a1, b1, c1, d1, cosine);
}

template < class K >
inline
typename K::Comparison_result
compare_dihedral_angle(const typename K::Vector_3& ab1, 
                       const typename K::Vector_3& ac1,
                       const typename K::Vector_3& ad1,
                       const typename K::Vector_3& ab2,
                       const typename K::Vector_3& ac2,
                       const typename K::Vector_3& ad2,
                       const K& k)
{
  return k.compare_dihedral_angle_3_object()(ab1, ac1, ad1, ab2, ac2, ad2);
}

template < class K >
inline
typename K::Comparison_result
compare_dihedral_angle(const typename K::Vector_3& ab1, 
                       const typename K::Vector_3& ac1,
                       const typename K::Vector_3& ad1,
                       const typename K::FT& cosine,
                       const K& k)
{
  return k.compare_dihedral_angle_3_object()(ab1, ac1, ad1, cosine);
}

template <class K, class T1, class T2, class T3>
inline
typename boost::enable_if<
  boost::mpl::equal_to<boost::mpl::integral_c<int,
                                              Ambient_dimension<T1>::type::value>,
                       boost::mpl::integral_c<int, 3> >,
  typename K::Comparison_result>
::type
  // boost::mpl::equal_to<typename Ambient_dimension<T1>::type,
  //                      boost::mpl::int_<3> >,
  // typename K::Comparison_result>::type
compare_distance(const T1 &o1,
                 const T2 &o2,
                 const T3 &o3, const K& k)
{
  return k.compare_distance_3_object()(o1, o2, o3);
}

template <class K, class T1, class T2, class T3, class T4>
inline
typename boost::enable_if<
  boost::mpl::equal_to<boost::mpl::integral_c<int,
                                              Ambient_dimension<T1>::type::value>,
                       boost::mpl::integral_c<int, 3> >,
  typename K::Comparison_result>
::type
  // boost::mpl::equal_to<typename Ambient_dimension<T1>::type,
  //                      boost::mpl::int_<3> >,
  // typename K::Comparison_result>::type
compare_distance(const T1 &o1,
                 const T2 &o2,
                 const T3 &o3,
                 const T4 &o4, const K& k)
{
  return k.compare_distance_3_object()(o1, o2, o3, o4);
}

template < class K >
inline
typename K::Comparison_result
compare_distance_to_point(const typename K::Point_3 &p,
                          const typename K::Point_3 &q,
                          const typename K::Point_3 &r,
			  const K& k)
{
  return k.compare_distance_3_object()(p, q, r);
}

template < class K >
inline
typename K::Comparison_result
compare_power_distance(const typename K::Point_3 &r,
                       const typename K::Weighted_point_3 &p,
                       const typename K::Weighted_point_3 &q,
                       const K& k)
{
  return k.compare_power_distance_3_object()(r, p, q);
}

template < class K >
inline
typename K::Comparison_result
compare_slope(const typename K::Point_3 &p,
              const typename K::Point_3 &q,
              const typename K::Point_3 &r,
              const typename K::Point_3 &s,
              const K& k)
{
  return k.compare_slope_3_object()(p, q, r, s);
}

template < class K >
inline
typename K::Comparison_result
compare_squared_distance(const typename K::Point_3 &p,
                         const typename K::Point_3 &q,
                         const typename K::FT &d2,
		         const K& k)
{
  return k.compare_squared_distance_3_object()(p, q, d2);
}

template < class K >
inline
typename K::Comparison_result
compare_squared_radius(const typename K::Point_3 &p,
		       const typename K::FT &sr,
		       const K& k)
{
  return k.compare_squared_radius_3_object()(p, sr);
}

template < class K >
inline
typename K::Comparison_result
compare_squared_radius(const typename K::Point_3 &p,
		       const typename K::Point_3 &q,
		       const typename K::FT &sr,
		       const K& k)
{
  return k.compare_squared_radius_3_object()(p, q, sr);
}

template < class K >
inline
typename K::Comparison_result
compare_squared_radius(const typename K::Point_3 &p,
		       const typename K::Point_3 &q,
		       const typename K::Point_3 &r,
		       const typename K::FT &sr,
		       const K& k)
{
  return k.compare_squared_radius_3_object()(p, q, r, sr);
}

template < class K >
inline
typename K::Comparison_result
compare_squared_radius(const typename K::Point_3 &p,
		       const typename K::Point_3 &q,
		       const typename K::Point_3 &r,
		       const typename K::Point_3 &s,
		       const typename K::FT &sr,
		       const K& k)
{
  return k.compare_squared_radius_3_object()(p, q, r, s, sr);
}

template < class K >
inline
typename K::Comparison_result
compare_lexicographically_xyz(const typename K::Point_3 &p,
                              const typename K::Point_3 &q,
			      const K& k)
{
  return k.compare_xyz_3_object()(p, q);
}

template < class K >
inline
typename K::Comparison_result
compare_signed_distance_to_plane(const typename K::Plane_3 &h,
				 const typename K::Point_3 &p,
				 const typename K::Point_3 &q,
				 const K &k)
{ 
  if (k.less_signed_distance_to_plane_3_object()(h, p, q)) return SMALLER;
  if (k.less_signed_distance_to_plane_3_object()(h, q, p)) return LARGER;
  return EQUAL;
}

template < class K >
inline
typename K::Comparison_result
compare_signed_distance_to_plane(const typename K::Point_3 &hp,
				 const typename K::Point_3 &hq,
				 const typename K::Point_3 &hr,
				 const typename K::Point_3 &p,
				 const typename K::Point_3 &q,
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
typename K::Comparison_result
compare_weighted_squared_radius(const typename K::Weighted_point_3 &p,
                                const typename K::FT &w, const K &k)
{
  return k.compare_weighted_squared_radius_3_object()(p, w);
}

template < class K >
inline
typename K::Comparison_result
compare_weighted_squared_radius(const typename K::Weighted_point_3 &p,
                                const typename K::Weighted_point_3 &q,
                                const typename K::FT &w, const K &k)
{
  return k.compare_weighted_squared_radius_3_object()(p, q, w);
}

template < class K >
inline
typename K::Comparison_result
compare_weighted_squared_radius(const typename K::Weighted_point_3 &p,
                                const typename K::Weighted_point_3 &q,
                                const typename K::Weighted_point_3 &r,
                                const typename K::FT &w, const K &k)
{
  return k.compare_weighted_squared_radius_3_object()(p, q, r, w);
}

template < class K >
inline
typename K::Comparison_result
compare_weighted_squared_radius(const typename K::Weighted_point_3 &p,
                                const typename K::Weighted_point_3 &q,
                                const typename K::Weighted_point_3 &r,
                                const typename K::Weighted_point_3 &s,
                                const typename K::FT &w, const K &k)
{
  return k.compare_weighted_squared_radius_3_object()(p, q, r, s, w);
}

template < class K >
inline
typename K::Comparison_result
compare_x(const typename K::Point_3 &p,
	  const typename K::Point_3 &q,
	  const K &k)
{ 
  return k.compare_x_3_object()(p, q);
}

template < class K >
inline
typename K::Comparison_result
compare_y(const typename K::Point_3 &p,
	  const typename K::Point_3 &q,
	  const K &k)
{ 
  return k.compare_y_3_object()(p, q);
}

template < class K >
inline
typename K::Comparison_result
compare_z(const typename K::Point_3 &p,
	  const typename K::Point_3 &q,
	  const K &k)
{ 
  return k.compare_z_3_object()(p, q);
}

template < class K >
inline
typename K::Comparison_result
compare_xyz(const typename K::Point_3 &p,
	    const typename K::Point_3 &q,
	    const K &k)
{ 
  return k.compare_xyz_3_object()(p, q);
}

template < class K >
inline
typename K::Boolean
coplanar(const typename K::Point_3 &p,
         const typename K::Point_3 &q,
         const typename K::Point_3 &r,
         const typename K::Point_3 &s, const K& k)
{
  return k.coplanar_3_object()(p, q, r, s);
}

template < class K >
inline
typename K::Orientation
coplanar_orientation(const typename K::Point_3 &p,
                     const typename K::Point_3 &q,
                     const typename K::Point_3 &r,
                     const typename K::Point_3 &s, const K& k)
{
  return k.coplanar_orientation_3_object()(p, q, r, s);
}

template < class K >
inline
typename K::Orientation
coplanar_orientation(const typename K::Point_3 &p,
                     const typename K::Point_3 &q,
                     const typename K::Point_3 &r, const K& k)
{
  return k.coplanar_orientation_3_object()(p, q, r);
}

template < class K >
inline
typename K::Bounded_side
coplanar_side_of_bounded_circle(const typename K::Point_3 &p,
                                const typename K::Point_3 &q,
                                const typename K::Point_3 &r,
                                const typename K::Point_3 &t,
                                const K& k)
{
  return k.coplanar_side_of_bounded_circle_3_object()(p, q, r, t);
}

template < class K >
inline
typename K::Vector_3
cross_product(const typename K::Vector_3 &v,
              const typename K::Vector_3 &w, const K& k)
{
  return k.construct_cross_product_vector_3_object()(v, w);
}

template < class K >
inline
typename K::FT
determinant(const typename K::Vector_3 &v0,
            const typename K::Vector_3 &v1,
            const typename K::Vector_3 &v2, const K &k)
{
  return k.compute_determinant_3_object()(v0, v1, v2);
}

template < class K >
inline
typename K::Line_3
equidistant_line(const typename K::Point_3 &p,
                 const typename K::Point_3 &q,
                 const typename K::Point_3 &r, const K& k)
{
  return k.construct_equidistant_line_3_object()(p, q, r);
}

template < class K >
inline
typename K::Boolean
has_smaller_distance_to_point(const typename K::Point_3 &p,
                              const typename K::Point_3 &q,
                              const typename K::Point_3 &r,
			      const K &k)
{
  return k.less_distance_to_point_3_object()(p, q, r);
}

template < class K >
inline
typename K::Boolean
has_larger_distance_to_point(const typename K::Point_3 &p,
			     const typename K::Point_3 &q,
			     const typename K::Point_3 &r,
			     const K &k)
{
  return k.compare_distance_3_object()(p, q, r) == LARGER;
}

template < class K >
inline
typename K::Boolean
has_larger_signed_distance_to_plane(const typename K::Plane_3 &h,
				    const typename K::Point_3 &p,
				    const typename K::Point_3 &q,
				    const K &k)
{ 
  return k.less_signed_distance_to_plane_3_object()(h, q, p);
}

template < class K >
inline
typename K::Boolean
has_larger_signed_distance_to_plane(const typename K::Point_3 &hp,
				    const typename K::Point_3 &hq,
				    const typename K::Point_3 &hr,
				    const typename K::Point_3 &p,
				    const typename K::Point_3 &q,
				    const K &k)
{ 
  return k.less_signed_distance_to_plane_3_object()(hp, hq, hr, q, p);
}

template < class K >
inline
typename K::Boolean
has_smaller_signed_distance_to_plane(const typename K::Plane_3 &h,
                                     const typename K::Point_3 &p,
                                     const typename K::Point_3 &q,
				     const K &k)
{ 
  return k.less_signed_distance_to_plane_3_object()(h, p, q);
}

template < class K >
inline
typename K::Boolean
has_smaller_signed_distance_to_plane(const typename K::Point_3 &hp,
                                     const typename K::Point_3 &hq,
                                     const typename K::Point_3 &hr,
                                     const typename K::Point_3 &p,
                                     const typename K::Point_3 &q,
				     const K &k)
{ 
  return k.less_signed_distance_to_plane_3_object()(hp, hq, hr, p, q);
}

template < class K >
inline
typename K::Boolean
less_x(const typename K::Point_3 &p,
       const typename K::Point_3 &q,
       const K &k)
{ 
  return k.less_x_3_object()(p, q);
}

template < class K >
inline
typename K::Boolean
less_y(const typename K::Point_3 &p,
       const typename K::Point_3 &q,
       const K &k)
{ 
  return k.less_y_3_object()(p, q);
}

template < class K >
inline
typename K::Boolean
less_z(const typename K::Point_3 &p,
       const typename K::Point_3 &q,
       const K &k)
{ 
  return k.less_z_3_object()(p, q);
}

template < class K >
inline
typename K::Boolean
lexicographically_xyz_smaller(const typename K::Point_3 &p,
                              const typename K::Point_3 &q,
                              const K &k)
{
  return k.less_xyz_3_object()(p, q);
}

template < class K >
inline
typename K::FT
l_infinity_distance(const typename K::Point_3 &p,
                    const typename K::Point_3 &q,
                    const K& k)
{
  return k.compute_L_infinity_distance_3_object()(p, q);
}

template < class K >
inline
typename K::Point_3
midpoint(const typename K::Point_3 &p,
         const typename K::Point_3 &q, const K &k)
{
  return k.construct_midpoint_3_object()(p, q);
}

template < class K >
inline
typename K::Point_3
max_vertex(const typename K::Iso_cuboid_3 &ic, const K &k)
{
  return k.construct_max_vertex_3_object()(ic);
}

template < class K >
inline
typename K::Point_3
min_vertex(const typename K::Iso_cuboid_3 &ic, const K &k)
{
  return k.construct_min_vertex_3_object()(ic);
}

template < class K >
inline
typename K::Vector_3
normal(const typename K::Point_3 &p, const typename K::Point_3 &q, const typename K::Point_3 &r, const K &k)
{
  return k.construct_normal_3_object()(p, q, r);
}
template < class K >
inline
typename K::Vector_3
unit_normal(const typename K::Point_3 &p, const typename K::Point_3 &q, const typename K::Point_3 &r, const K &k)
{
  return k.construct_unit_normal_3_object()(p, q, r);
}

template <class K >
inline
typename K::Orientation
orientation(const typename K::Point_3 &p,
	    const typename K::Point_3 &q,
	    const typename K::Point_3 &r,
	    const typename K::Point_3 &s, const K &k)
{
  return k.orientation_3_object()(p, q, r, s);
}

template <class K >
inline
typename K::Orientation
orientation(const typename K::Vector_3 &u,
	    const typename K::Vector_3 &v,
	    const typename K::Vector_3 &w, const K &k)
{
  return k.orientation_3_object()(u, v, w);
}

template < class K >
inline
typename K::Vector_3
orthogonal_vector(const typename K::Point_3 &p,
		  const typename K::Point_3 &q,
		  const typename K::Point_3 &r, const K &k)
{
  return k.construct_orthogonal_vector_3_object()(p, q, r);
}

template < class K >
inline
typename K::Vector_3
orthogonal_vector(const typename K::Plane_3 &p, const K &k)
{
  return k.construct_orthogonal_vector_3_object()(p);
}

template <typename K>
inline
typename K::Boolean
parallel(const typename K::Line_3 &l1,
         const typename K::Line_3 &l2, const K &k)
{
  return k.are_parallel_3_object()(l1, l2);
}

template <typename K>
inline
typename K::Boolean
parallel(const typename K::Plane_3 &h1,
         const typename K::Plane_3 &h2, const K &k)
{
  return k.are_parallel_3_object()(h1, h2);
}

template <typename K>
inline
typename K::Boolean
parallel(const typename K::Ray_3 &r1,
         const typename K::Ray_3 &r2, const K &k)
{
  return k.are_parallel_3_object()(r1, r2);
}

template <typename K>
inline
typename K::Boolean
parallel(const typename K::Segment_3 &s1,
         const typename K::Segment_3 &s2, const K &k)
{
  return k.are_parallel_3_object()(s1, s2);
}

template <typename K>
inline
typename K::FT
power_distance_to_power_sphere(const typename K::Weighted_point_3 &p,
                               const typename K::Weighted_point_3 &q,
                               const typename K::Weighted_point_3 &r,
                               const typename K::Weighted_point_3 &s,
                               const typename K::Weighted_point_3 &t, const K &k)
{
  return k.compute_power_distance_to_power_sphere_3_object()(p, q, r, s, t);
}

template <class K >
inline
typename K::FT
power_product(const typename K::Weighted_point_3 &p,
              const typename K::Weighted_point_3 &q, const K &k)
{
  return k.compute_power_product_3_object()(p, q);
}

template <class K >
inline
typename K::Bounded_side
power_side_of_bounded_power_sphere(const typename K::Weighted_point_3 &p,
                                   const typename K::Weighted_point_3 &q, const K &k)
{
  return k.power_side_of_bounded_power_sphere_3_object()(p, q);
}

template <class K >
inline
typename K::Bounded_side
power_side_of_bounded_power_sphere(const typename K::Weighted_point_3 &p,
                                   const typename K::Weighted_point_3 &q,
                                   const typename K::Weighted_point_3 &r, const K &k)
{
  return k.power_side_of_bounded_power_sphere_3_object()(p, q, r);
}

template <class K >
inline
typename K::Bounded_side
power_side_of_bounded_power_sphere(const typename K::Weighted_point_3 &p,
                                   const typename K::Weighted_point_3 &q,
                                   const typename K::Weighted_point_3 &r,
                                   const typename K::Weighted_point_3 &s, const K &k)
{
  return k.power_side_of_bounded_power_sphere_3_object()(p, q, r, s);
}

template <class K >
inline
typename K::Bounded_side
power_side_of_bounded_power_sphere(const typename K::Weighted_point_3 &p,
                                   const typename K::Weighted_point_3 &q,
                                   const typename K::Weighted_point_3 &r,
                                   const typename K::Weighted_point_3 &s,
                                   const typename K::Weighted_point_3 &t, const K &k)
{
  return k.power_side_of_bounded_power_sphere_3_object()(p, q, r, s, t);
}

template <class K >
inline
typename K::Oriented_side
power_side_of_oriented_power_sphere(const typename K::Weighted_point_3 &p,
                                    const typename K::Weighted_point_3 &q, const K &k)
{
  return k.power_side_of_oriented_power_sphere_3_object()(p, q);
}

template <class K >
inline
typename K::Oriented_side
power_side_of_oriented_power_sphere(const typename K::Weighted_point_3 &p,
                                    const typename K::Weighted_point_3 &q,
                                    const typename K::Weighted_point_3 &r, const K &k)
{
  return k.power_side_of_oriented_power_sphere_3_object()(p, q, r);
}

template <class K >
inline
typename K::Oriented_side
power_side_of_oriented_power_sphere(const typename K::Weighted_point_3 &p,
                                    const typename K::Weighted_point_3 &q,
                                    const typename K::Weighted_point_3 &r,
                                    const typename K::Weighted_point_3 &s, const K &k)
{
  return k.power_side_of_oriented_power_sphere_3_object()(p, q, r, s);
}

template <class K >
inline
typename K::Oriented_side
power_side_of_oriented_power_sphere(const typename K::Weighted_point_3 &p,
                                    const typename K::Weighted_point_3 &q,
                                    const typename K::Weighted_point_3 &r,
                                    const typename K::Weighted_point_3 &s,
                                    const typename K::Weighted_point_3 &t, const K &k)
{
  return k.power_side_of_oriented_power_sphere_3_object()(p, q, r, s, t);
}

template <class K >
inline
typename K::Bounded_side
side_of_bounded_sphere(const typename K::Point_3 &p,
                       const typename K::Point_3 &q,
                       const typename K::Point_3 &test, const K &k)
{
  return k.side_of_bounded_sphere_3_object()(p, q, test);
}

template <class K >
inline
typename K::Bounded_side
side_of_bounded_sphere(const typename K::Point_3 &p,
                       const typename K::Point_3 &q,
                       const typename K::Point_3 &r,
                       const typename K::Point_3 &test, const K &k)
{
  return k.side_of_bounded_sphere_3_object()(p, q, r, test);
}

template <class K >
inline
typename K::Bounded_side
side_of_bounded_sphere(const typename K::Point_3 &p,
                       const typename K::Point_3 &q,
                       const typename K::Point_3 &r,
                       const typename K::Point_3 &s,
                       const typename K::Point_3 &test, const K &k)
{
  return k.side_of_bounded_sphere_3_object()(p, q, r, s, test);
}

template <class K >
inline
typename K::Oriented_side
side_of_oriented_sphere(const typename K::Point_3 &p,
                        const typename K::Point_3 &q,
                        const typename K::Point_3 &r,
                        const typename K::Point_3 &s,
                        const typename K::Point_3 &test, const K &k)
{
  return k.side_of_oriented_sphere_3_object()(p, q, r, s, test);
}

template <typename K>
inline
typename K::FT
squared_area(const typename K::Point_3 &p,
	     const typename K::Point_3 &q,
	     const typename K::Point_3 &r, const K &k)
{
  return k.compute_squared_area_3_object()(p, q, r);
}

template < class K >
inline
typename K::FT
squared_radius(const typename K::Point_3 &p,
	       const typename K::Point_3 &q,
	       const typename K::Point_3 &r,
	       const typename K::Point_3 &s, const K &k)
{
  return k.compute_squared_radius_3_object()(p, q, r, s);
}

template < class K >
inline
typename K::FT
squared_radius(const typename K::Point_3 &p,
	       const typename K::Point_3 &q,
	       const typename K::Point_3 &r, const K &k)
{
  return k.compute_squared_radius_3_object()(p, q, r);
}

template < class K >
inline
typename K::FT
squared_radius(const typename K::Point_3 &p,
	       const typename K::Point_3 &q, const K &k)
{
  return k.compute_squared_radius_3_object()(p, q);
}

template < class K >
inline
typename K::FT
squared_radius(const typename K::Point_3 &p, const K &k)
{
  return k.compute_squared_radius_3_object()(p);
}

template < class K >
inline
typename K::FT
squared_radius_smallest_orthogonal_sphere(const typename K::Weighted_point_3 &p,
                                          const K &k)
{
  return k.compute_squared_radius_smallest_orthogonal_sphere_3_object()(p);
}

template < class K >
inline
typename K::FT
squared_radius_smallest_orthogonal_sphere(const typename K::Weighted_point_3 &p,
                                          const typename K::Weighted_point_3 &q,
                                          const K &k)
{
  return k.compute_squared_radius_smallest_orthogonal_sphere_3_object()(p, q);
}

template < class K >
inline
typename K::FT
squared_radius_smallest_orthogonal_sphere(const typename K::Weighted_point_3 &p,
                                          const typename K::Weighted_point_3 &q,
                                          const typename K::Weighted_point_3 &r,
                                          const K &k)
{
  return k.compute_squared_radius_smallest_orthogonal_sphere_3_object()(p, q, r);
}

template < class K >
inline
typename K::FT
squared_radius_smallest_orthogonal_sphere(const typename K::Weighted_point_3 &p,
                                          const typename K::Weighted_point_3 &q,
                                          const typename K::Weighted_point_3 &r,
                                          const typename K::Weighted_point_3 &s,
                                          const K &k)
{
  return k.compute_squared_radius_smallest_orthogonal_sphere_3_object()(p, q, r, s);
}

template < class K >
inline
typename K::FT
volume(const typename K::Point_3 &p,
       const typename K::Point_3 &q,
       const typename K::Point_3 &r,
       const typename K::Point_3 &s, const K &k)
{
  return k.compute_volume_3_object()(p, q, r, s);
}

template < class K >
inline
typename K::Point_3
weighted_circumcenter(const typename K::Weighted_point_3 &p,
                      const typename K::Weighted_point_3 &q, const K &k)
{
  return k.construct_weighted_circumcenter_3_object()(p, q);
}

template < class K >
inline
typename K::Point_3
weighted_circumcenter(const typename K::Weighted_point_3 &p,
                      const typename K::Weighted_point_3 &q,
                      const typename K::Weighted_point_3 &r, const K &k)
{
  return k.construct_weighted_circumcenter_3_object()(p, q, r);
}

template < class K >
inline
typename K::Point_3
weighted_circumcenter(const typename K::Weighted_point_3 &p,
                      const typename K::Weighted_point_3 &q,
                      const typename K::Weighted_point_3 &r,
                      const typename K::Weighted_point_3 &s, const K &k)
{
  return k.construct_weighted_circumcenter_3_object()(p, q, r, s);
}

template < class K >
inline
typename K::Boolean
x_equal(const typename K::Point_3 &p,
        const typename K::Point_3 &q, const K &k)
{
  return k.equal_x_3_object()(p, q);
}

template < class K >
inline
typename K::Boolean
y_equal(const typename K::Point_3 &p,
        const typename K::Point_3 &q, const K &k)
{
  return k.equal_y_3_object()(p, q);
}

template < class K >
inline
typename K::Boolean
z_equal(const typename K::Point_3 &p,
        const typename K::Point_3 &q, const K &k)
{
  return k.equal_z_3_object()(p, q);
}


// The following functions only call some of the previous ones.

template <typename K>
inline
typename K::Boolean
are_negative_oriented(const typename K::Point_3 &p,
                      const typename K::Point_3 &q,
                      const typename K::Point_3 &r,
                      const typename K::Point_3 &s, const K &k)
{
  return internal::orientation(p, q, r, s, k) == NEGATIVE;
}

template <typename K>
inline
typename K::Boolean
are_positive_oriented(const typename K::Point_3 &p,
                      const typename K::Point_3 &q,
                      const typename K::Point_3 &r,
                      const typename K::Point_3 &s, const K &k)
{
  return internal::orientation(p, q, r, s, k) == POSITIVE;
}

template < class K >
inline
typename K::Boolean
lexicographically_xyz_smaller_or_equal(const typename K::Point_3 &p,
                                       const typename K::Point_3 &q,
                                       const K&k)
{
  return internal::compare_lexicographically_xyz(p, q, k) != LARGER;
}

} // namespace internal

} //namespace CGAL

#endif  // CGAL_KERNEL_GLOBAL_FUNCTIONS_INTERNAL_3_H
