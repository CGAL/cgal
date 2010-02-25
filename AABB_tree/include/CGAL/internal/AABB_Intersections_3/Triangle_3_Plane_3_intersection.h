// Copyright (c) 2003  INRIA Sophia-Antipolis (France).
// Copyright (c) 2009  GeometryFactory (France).
// All rights reserved.
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
// Adapted from <CGAL/Triangle_3_Plane_3_do_intersect.h>
//
// Author(s)     : Philippe Guigue, Laurent Rineau

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_TRIANGLE_3_PLANE_3_INTERSECTION_H
#define CGAL_INTERNAL_INTERSECTIONS_3_TRIANGLE_3_PLANE_3_INTERSECTION_H

#include <CGAL/Plane_3.h>
#include <CGAL/Triangle_3.h>

namespace CGAL {

namespace internal {

template <class K>
inline
typename K::Point_3
inter_plane_triangle_3_aux(const typename K::Point_3 &p1,
                           const typename K::RT & f1,
                           const typename K::Point_3 &p2,
                           const typename K::RT & f2)
{
  return typename K::Point_3(f2 * p1.x() - f1 * p2.x(),
                             f2 * p1.y() - f1 * p2.y(),
                             f2 * p1.z() - f1 * p2.z(),
                             f2 - f1);
}

template <class K>
Object
intersection(const typename K::Plane_3  &plane,
	     const typename K::Triangle_3 &triangle,
	     const K& k)
{
  CGAL_kernel_precondition( ! k.is_degenerate_3_object()(triangle)) ;
  CGAL_kernel_precondition( ! k.is_degenerate_3_object()(plane)) ;

  typedef typename K::Point_3 Point_3;
  typedef typename K::Segment_3 Segment_3;
  typedef typename K::Object_3 Object_3;
  typedef typename K::RT RT;
  typedef typename Sgn<RT>::result_type SignRT;

  typename K::Construct_vertex_3 vertex_on =
    k.construct_vertex_3_object();

  //typename K::Oriented_side_3 oriented_side =
    //k.oriented_side_3_object(); // PA: never used

  typename K::Construct_segment_3 segment =
    k.construct_segment_3_object();

  typename K::Construct_object_3 make_object =
    k.construct_object_3_object();

  const Point_3& t0 = vertex_on(triangle,0);
  const Point_3& t1 = vertex_on(triangle,1);
  const Point_3& t2 = vertex_on(triangle,2);

  const RT f0 = plane.a()*t0.hx() + plane.b()*t0.hy()
    + plane.c()*t0.hz() + wmult_hw((K*)0, plane.d(), t0);
  const RT f1 = plane.a()*t1.hx() + plane.b()*t1.hy()
    + plane.c()*t1.hz() + wmult_hw((K*)1, plane.d(), t1);
  const RT f2 = plane.a()*t2.hx() + plane.b()*t2.hy()
    + plane.c()*t2.hz() + wmult_hw((K*)2, plane.d(), t2);

  const SignRT s0 = CGAL_NTS sign(f0);
  const SignRT s1 = CGAL_NTS sign(f1);
  const SignRT s2 = CGAL_NTS sign(f2);

  if(s0 == ZERO) {
    if(s1 == ZERO) {
      if(s2 == ZERO) {
        // all zero
        return make_object(triangle);
      }
      else {
        // s0, s1 zero
        return make_object(segment(t0, t1));
      }
    }
    else if(s2 == ZERO) {
      // s0, s2 zero
      return make_object(segment(t0, t2));
    }
    else {
      // s0 zero
      if(s1 == s2){
        return make_object(t0);
      } else {
        return make_object(segment(t0,
                                   inter_plane_triangle_3_aux<K>(t1, f1,
                                                                 t2, f2)));
      }
    }
  }
  else {
    if(s1 == ZERO) {
      if(s2 == ZERO) {
        // s1, s2 zero
        return make_object(segment(t1, t2));
      }
      else {
        // s1 zero
        if(s0 == s2){
          return make_object(t1);
        } else {
          return make_object(segment(t1,
                                     inter_plane_triangle_3_aux<K>(t0, f0,
                                                                   t2, f2)));
        }
      }
    }
    else if( s2 == ZERO ) {
      // s2 zero
      if(s0 == s1){
        return make_object(t2);
      } else {
        return make_object(segment(t2,
                                   inter_plane_triangle_3_aux<K>(t0, f0,
                                                              t1, f1)));
      }
    }
    else {
      // all non-zero
      if(s0 == s1) {
        if(s0 == s2) {
          return Object_3();
        }
        else {
          // s0 and s1 on same side, s2 different
          return make_object(segment(inter_plane_triangle_3_aux<K>(t0, f0,
                                                                t2, f2),
                                     inter_plane_triangle_3_aux<K>(t1, f1,
                                                                t2, f2)));
        }
      }
      else if(s0 == s2) {
        // s0 and s2 on same side, s1 different
          return make_object(segment(inter_plane_triangle_3_aux<K>(t1, f1,
                                                                t2, f2),
                                     inter_plane_triangle_3_aux<K>(t1, f1,
                                                                t0, f0)));
      }
      else {
        // s1 and s2 on same side, s0 different
          return make_object(segment(inter_plane_triangle_3_aux<K>(t0, f0,
                                                                t2, f2),
                                     inter_plane_triangle_3_aux<K>(t0, f0,
                                                                t1, f1)));
      }
    }
  }
} // end function internal::intersection(Plane, Triangle)

template <class K>
inline
Object
intersection(const typename K::Triangle_3 &triangle,
	     const typename K::Plane_3  &plane,
	     const K& k)
{
  return intersection(plane, triangle, k);
}

} // end namespace internal

template <class K>
inline
Object
intersection(const Plane_3<K> &plane, const Triangle_3<K> &triangle)
{
  return typename K::Intersect_3()(plane, triangle);
}

template <class K>
inline
Object
intersection(const Triangle_3<K> &triangle, const Plane_3<K> &plane)
{
  return typename K::Intersect_3()(plane, triangle);
}

} // end namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_TRIANGLE_3_PLANE_3_INTERSECTION_H
