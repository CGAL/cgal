// Copyright (c) 2009  GeometryFactory (France), INRIA Sophia-Antipolis (France).
// All rights reserved.
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
// Author(s)     :  Laurent Rineau, Stephane Tayeb
//
// Note: This implementation is adapted from Triangle_3_Segment_3_do_intersect.h

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_TRIANGLE_3_SEGMENT_3_INTERSECTION_H
#define CGAL_INTERNAL_INTERSECTIONS_3_TRIANGLE_3_SEGMENT_3_INTERSECTION_H

#include <CGAL/kernel_basic.h>
#include <CGAL/intersections.h>

namespace CGAL {
namespace internal {

template <class K>
typename K::Point_3
t3s3_intersection_coplanar_aux(const typename K::Point_3& p,
                               const typename K::Point_3& q,
                               const typename K::Point_3& a,
                               const typename K::Point_3& b,
                               const K& k)
{
  // Returns the intersection point between segment [p,q] and [a,b]
  //
  // preconditions:
  //   + p,q,a,b are coplanar

  typedef typename K::Vector_3 Vector_3;
  typedef typename K::FT FT;

  typename K::Construct_vector_3 vector =
    k.construct_vector_3_object();

  typename K::Construct_cross_product_vector_3 cross_product =
    k.construct_cross_product_vector_3_object();

  typename K::Compute_scalar_product_3 scalar_product =
    k.compute_scalar_product_3_object();

  typename K::Compute_squared_length_3 sq_length =
    k.compute_squared_length_3_object();

  const Vector_3 pq = vector(p,q);
  const Vector_3 ab = vector(a,b);
  const Vector_3 pa = vector(p,a);

  const Vector_3 pa_ab = cross_product(pa,ab);
  const Vector_3 pq_ab = cross_product(pq,ab);

  const FT t = scalar_product(pa_ab,pq_ab) / sq_length(pq_ab);

  return ( p + t*pq );
}


template <class K>
typename Intersection_traits<K, typename K::Triangle_3, typename K::Segment_3>::result_type
t3s3_intersection_coplanar_aux(const typename K::Point_3& a,
                               const typename K::Point_3& b,
                               const typename K::Point_3& c,
                               const typename K::Point_3& p,
                               const typename K::Point_3& q,
                               const bool negative_side,
                               const K& k)
{
  // This function is designed to clip pq into the triangle abc.
  // Point configuration should be as follows
  //
  //     +p
  //     |    +b
  //     |
  //  +c |       +a
  //     |
  //     +q
  //
  // We know that c is isolated on the negative side of pq, but we don't know
  // p position wrt [bc]&[ca] and q position wrt [bc]&[ca]

  typedef typename K::Point_3 Point_3;

  typename K::Coplanar_orientation_3 coplanar_orientation =
    k.coplanar_orientation_3_object();

  typename K::Construct_segment_3 segment =
    k.construct_segment_3_object();

  const Orientation bcq = coplanar_orientation(b,c,q);
  const Orientation cap = coplanar_orientation(c,a,p);

  if ( NEGATIVE == bcq || NEGATIVE == cap )
    return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Segment_3>();
  else if ( COLLINEAR == bcq )
    // q is inside [c,b], p is outside t (because of pqc)
    return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Segment_3>(q);
  else if ( COLLINEAR == cap )
    // p is inside [c,a], q is outside t (because of pqc)
    return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Segment_3>(p);
  else // bcq == POSITIVE && cap == POSITIVE
  {
    // Here we know the intersection is not empty
    // Let's get the intersection points
    Point_3 p_side_end_point(p);
    if ( NEGATIVE == coplanar_orientation(b,c,p) )
    {
      p_side_end_point = t3s3_intersection_coplanar_aux(p,q,b,c,k);
    }

    Point_3 q_side_end_point(q);
    if ( NEGATIVE == coplanar_orientation(c,a,q) )
    {
      q_side_end_point = t3s3_intersection_coplanar_aux(p,q,c,a,k);
    }

    if ( negative_side )
      return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Segment_3>(segment(p_side_end_point, q_side_end_point));
    else
      return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Segment_3>(segment(q_side_end_point, p_side_end_point));
  }
}


template <class K>
typename Intersection_traits<K, typename K::Triangle_3, typename K::Segment_3>::result_type
t3s3_intersection_collinear_aux(const typename K::Point_3& a,
                               const typename K::Point_3& b,
                               const typename K::Point_3& p,
                               const typename K::Point_3& q,
                               const K& k)
{
  // Builds resulting segment of intersection of [a,b] and [p,q]
  // Precondition: [a,b] and [p,q] have the same direction
  typename K::Construct_segment_3 segment =
    k.construct_segment_3_object();

  typename K::Collinear_are_ordered_along_line_3 collinear_ordered =
    k.collinear_are_ordered_along_line_3_object();

  typename K::Equal_3 equals = k.equal_3_object();
 
  // possible orders: [p,a,b,q], [p,a,q,b], [p,q,a,b], [a,p,b,q], [a,p,q,b], [a,b,p,q]
  if ( collinear_ordered(p,a,b) )
  {
    // p is before a
    //possible orders: [p,a,b,q], [p,a,q,b], [p,q,a,b]
    if ( collinear_ordered(a,b,q) )
      return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Segment_3>(segment(a,b)); //[p,a,b,q]
    else{
      if ( collinear_ordered(q,a,b) )
        return equals(a,q)? //[p,q,a,b]
             intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Segment_3>(a):
             intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Segment_3>();
      return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Segment_3>(segment(a,q));
    }
  }
  else
  {
    // p is after a
    //possible orders: [a,p,b,q], [a,p,q,b], [a,b,p,q]
    if ( collinear_ordered(p,b,q) )
      return equals(p,b)? // [a,p,b,q]
             intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Segment_3>(p):
             intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Segment_3>(segment(p,b));
    else{
      if ( collinear_ordered(a,b,p) )
        return equals(p,b)? // [a,b,p,q]
          intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Segment_3>(p) :
          intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Segment_3>();
      return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Segment_3>(segment(p,q)); // [a,p,q,b]
    }
  }
}


template <class K>
typename Intersection_traits<K, typename K::Segment_3, typename K::Triangle_3>::result_type
intersection_coplanar(const typename K::Triangle_3 &t,
                      const typename K::Segment_3  &s,
                      const K & k )
{

  CGAL_kernel_precondition( ! k.is_degenerate_3_object()(t) ) ;
  CGAL_kernel_precondition( ! k.is_degenerate_3_object()(s) ) ;

  typedef typename K::Point_3 Point_3;

  typename K::Construct_point_on_3 point_on =
    k.construct_point_on_3_object();

  typename K::Construct_vertex_3 vertex_on =
    k.construct_vertex_3_object();

  typename K::Coplanar_orientation_3 coplanar_orientation =
    k.coplanar_orientation_3_object();

  typename K::Collinear_are_ordered_along_line_3 collinear_ordered =
    k.collinear_are_ordered_along_line_3_object();

  const Point_3 & p = point_on(s,0);
  const Point_3 & q = point_on(s,1);

  const Point_3 & A = vertex_on(t,0);
  const Point_3 & B = vertex_on(t,1);
  const Point_3 & C = vertex_on(t,2);

  int k0 = 0;
  int k1 = 1;
  int k2 = 2;

  // Determine the orientation of the triangle in the common plane
  if (coplanar_orientation(A,B,C) != POSITIVE)
  {
    // The triangle is not counterclockwise oriented
    // swap two vertices.
    std::swap(k1,k2);
  }

  const Point_3& a = vertex_on(t,k0);
  const Point_3& b = vertex_on(t,k1);
  const Point_3& c = vertex_on(t,k2);

  // Test whether the segment's supporting line intersects the
  // triangle in the common plane
  const Orientation pqa = coplanar_orientation(p,q,a);
  const Orientation pqb = coplanar_orientation(p,q,b);
  const Orientation pqc = coplanar_orientation(p,q,c);

  switch ( pqa ) {
    // -----------------------------------
    // pqa POSITIVE
    // -----------------------------------
    case POSITIVE:
      switch ( pqb ) {
        case POSITIVE:
          switch ( pqc ) {
            case POSITIVE:
              // the triangle lies in the positive halfspace
              // defined by the segment's supporting line.
              return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Triangle_3>();
            case NEGATIVE:
              // c is isolated on the negative side
              return t3s3_intersection_coplanar_aux(a,b,c,p,q,true,k);
            default: // COLLINEAR
              if ( collinear_ordered(p,c,q) ) // c is inside [p,q]
                return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Triangle_3>(c);
              else
                return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Triangle_3>();
          }
        case NEGATIVE:
          if ( POSITIVE == pqc )
            // b is isolated on the negative side
            return t3s3_intersection_coplanar_aux(c,a,b,p,q,true,k);
          else
            // a is isolated on the positive side
            return t3s3_intersection_coplanar_aux(b,c,a,q,p,false,k);
        case COLLINEAR:
          switch ( pqc ) {
            case POSITIVE:
              if ( collinear_ordered(p,b,q) ) // b is inside [p,q]
                return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Triangle_3>(b);
              else
                return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Triangle_3>();
            case NEGATIVE:
              // a is isolated on the positive side
              return t3s3_intersection_coplanar_aux(b,c,a,q,p,false,k);
            default: // COLLINEAR
              // b,c,p,q are aligned, [p,q]&[b,c] have the same direction
              return t3s3_intersection_collinear_aux(b,c,p,q,k);
          }
        default: // should not happen.
          CGAL_error();
          return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Triangle_3>();
      }

    // -----------------------------------
    // pqa NEGATIVE
    // -----------------------------------
    case NEGATIVE:
      switch ( pqb ) {
        case POSITIVE:
          if ( POSITIVE == pqc )
            // a is isolated on the negative side
            return t3s3_intersection_coplanar_aux(b,c,a,p,q,true,k);
          else
            // b is isolated on the positive side
            return t3s3_intersection_coplanar_aux(c,a,b,q,p,false,k);

        case NEGATIVE:
          switch ( pqc ) {
            case POSITIVE:
              // c is isolated on the positive side
              return t3s3_intersection_coplanar_aux(a,b,c,q,p,false,k);
            case NEGATIVE:
              // the triangle lies in the negative halfspace
              // defined by the segment's supporting line.
              return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Triangle_3>();
            default: // COLLINEAR
              if ( collinear_ordered(p,c,q) ) // c is inside [p,q]
                return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Triangle_3>(c);
              else
                return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Triangle_3>();
          }

        case COLLINEAR:
          switch ( pqc ) {
            case POSITIVE:
              // a is isolated on the negative side
              return t3s3_intersection_coplanar_aux(b,c,a,p,q,true,k);
            case NEGATIVE:
              if ( collinear_ordered(p,b,q) ) // b is inside [p,q]
                return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Triangle_3>(b);
              else
                return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Triangle_3>();
            default: // COLLINEAR
              // b,c,p,q are aligned, [p,q]&[c,b] have the same direction
              return t3s3_intersection_collinear_aux(c,b,p,q,k);
          }

        default: // should not happen.
          CGAL_error();
          return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Triangle_3>();
      }

    // -----------------------------------
    // pqa COLLINEAR
    // -----------------------------------
    case COLLINEAR:
      switch ( pqb ) {
        case POSITIVE:
          switch ( pqc ) {
            case POSITIVE:
              if ( collinear_ordered(p,a,q) ) // a is inside [p,q]
                return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Triangle_3>(a);
              else
                return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Triangle_3>();
            case NEGATIVE:
              // b is isolated on the positive side
              return t3s3_intersection_coplanar_aux(c,a,b,q,p,false,k);
            default: // COLLINEAR
              // a,c,p,q are aligned, [p,q]&[c,a] have the same direction
              return t3s3_intersection_collinear_aux(c,a,p,q,k);
          }

        case NEGATIVE:
          switch ( pqc ) {
            case POSITIVE:
              // b is isolated on the negative side
              return t3s3_intersection_coplanar_aux(c,a,b,p,q,true,k);
            case NEGATIVE:
              if ( collinear_ordered(p,a,q) ) // a is inside [p,q]
                return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Triangle_3>(a);
              else
                return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Triangle_3>();
            default: // COLLINEAR
              // a,c,p,q are aligned, [p,q]&[a,c] have the same direction
              return t3s3_intersection_collinear_aux(a,c,p,q,k);
          }

        case COLLINEAR:
          switch ( pqc ) {
            case POSITIVE:
              // a,b,p,q are aligned, [p,q]&[a,b] have the same direction
              return t3s3_intersection_collinear_aux(a,b,p,q,k);
            case NEGATIVE:
              // a,b,p,q are aligned, [p,q]&[b,a] have the same direction
              return t3s3_intersection_collinear_aux(b,a,p,q,k);
            default:  // COLLINEAR
              // case pqc == COLLINEAR is impossible since the triangle is
              // assumed to be non flat
              CGAL_error();
              return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Triangle_3>();
          }

        default: // should not happen.
          CGAL_error();
          return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Triangle_3>();
      }

    default:// should not happen.
      CGAL_error();
      return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Triangle_3>();
  }
}


template <class K>
typename Intersection_traits<K, typename K::Segment_3, typename K::Triangle_3>::result_type
intersection(const typename K::Triangle_3 &t,
             const typename K::Segment_3  &s,
             const K & k)
{
  CGAL_kernel_precondition( ! k.is_degenerate_3_object()(t) ) ;
  CGAL_kernel_precondition( ! k.is_degenerate_3_object()(s) ) ;

  typedef typename K::Point_3 Point_3;

  typename K::Construct_point_on_3 point_on =
    k.construct_point_on_3_object();

  typename K::Construct_vertex_3 vertex_on =
    k.construct_vertex_3_object();

  typename K::Orientation_3 orientation =
    k.orientation_3_object();

  const Point_3 & a = vertex_on(t,0);
  const Point_3 & b = vertex_on(t,1);
  const Point_3 & c = vertex_on(t,2);
  const Point_3 & p = point_on(s,0);
  const Point_3 & q = point_on(s,1);

  const Orientation abcp = orientation(a,b,c,p);
  const Orientation abcq = orientation(a,b,c,q);

  switch ( abcp ) {
    case POSITIVE:
      switch ( abcq ) {
        case POSITIVE:
          // the segment lies in the positive open halfspaces defined by the
          // triangle's supporting plane
          return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Triangle_3>();

        case NEGATIVE:
          // p sees the triangle in counterclockwise order
          if ( orientation(p,q,a,b) != POSITIVE
              && orientation(p,q,b,c) != POSITIVE
              && orientation(p,q,c,a) != POSITIVE )
          {
            // The intersection should be a point
            
            typename Intersection_traits<K, typename K::Line_3, typename K::Plane_3>
              ::result_type 
              v = internal::intersection(s.supporting_line(),t.supporting_plane(), K());
            if(v) {
              if(const Point_3* res = intersect_get<Point_3>(v))
                return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Triangle_3>(*res);
            else
                return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Triangle_3>();
            } else {
              return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Triangle_3>();
            }
          }
          else
            return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Triangle_3>();

        case COPLANAR:
          // q belongs to the triangle's supporting plane
          // p sees the triangle in counterclockwise order
          if ( orientation(p,q,a,b) != POSITIVE
              && orientation(p,q,b,c) != POSITIVE
              && orientation(p,q,c,a) != POSITIVE )
            return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Triangle_3>(q);
          else
            return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Triangle_3>();

        default: // should not happen.
          CGAL_error();
          return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Triangle_3>();
      }
    case NEGATIVE:
      switch ( abcq ) {
        case POSITIVE:
          // q sees the triangle in counterclockwise order
          if ( orientation(q,p,a,b) != POSITIVE
            && orientation(q,p,b,c) != POSITIVE
            && orientation(q,p,c,a) != POSITIVE )
          {
            typename Intersection_traits<K, typename K::Line_3, typename K::Plane_3>::result_type 
              v = internal::intersection(s.supporting_line(),t.supporting_plane(), K());
            if(v) {
              if(const Point_3* res = intersect_get<Point_3>(v))
                return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Triangle_3>(*res);
            else
                return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Triangle_3>();
            } else {
              return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Triangle_3>();
          }
          } else
            return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Triangle_3>();
        case NEGATIVE:
          // the segment lies in the negative open halfspaces defined by the
          // triangle's supporting plane
          return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Triangle_3>();

        case COPLANAR:
          // q belongs to the triangle's supporting plane
          // p sees the triangle in clockwise order
          if ( orientation(q,p,a,b) != POSITIVE
              && orientation(q,p,b,c) != POSITIVE
              && orientation(q,p,c,a) != POSITIVE )
            return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Triangle_3>(q);
          else
            return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Triangle_3>();

        default: // should not happen.
          CGAL_error();
          return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Triangle_3>();
      }
    case COPLANAR: // p belongs to the triangle's supporting plane
      switch ( abcq ) {
        case POSITIVE:
          // q sees the triangle in counterclockwise order
          if ( orientation(q,p,a,b) != POSITIVE
              && orientation(q,p,b,c) != POSITIVE
              && orientation(q,p,c,a) != POSITIVE )
            return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Triangle_3>(p);
          else
            return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Triangle_3>();

        case NEGATIVE:
          // q sees the triangle in clockwise order
          if ( orientation(p,q,a,b) != POSITIVE
              && orientation(p,q,b,c) != POSITIVE
              && orientation(p,q,c,a) != POSITIVE )
            return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Triangle_3>(p);
          else
            return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Triangle_3>();

        case COPLANAR:
          // the segment is coplanar with the triangle's supporting plane
          // we test whether the segment intersects the triangle in the common
          // supporting plane
          return intersection_coplanar(t,s,k);

        default: // should not happen.
          CGAL_error();
          return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Triangle_3>();
      }
    default: // should not happen.
      CGAL_error();
      return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Triangle_3>();
  }
}

template <class K>
inline
typename Intersection_traits<K, typename K::Segment_3, 
                               typename K::Triangle_3>::result_type
intersection(const typename K::Segment_3  &s,
             const typename K::Triangle_3 &t,
             const K & k)
{
  return internal::intersection(t,s,k);
}


} // end namespace internal

CGAL_INTERSECTION_FUNCTION(Triangle_3, Segment_3, 3)

} // end namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_TRIANGLE_3_SEGMENT_3_INTERSECTION_H
