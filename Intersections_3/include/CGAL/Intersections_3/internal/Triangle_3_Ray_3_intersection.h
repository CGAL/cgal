// Copyright (c) 2009  GeometryFactory (France), INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     :  Laurent Rineau, Stephane Tayeb
//
// Note: This implementation is adapted from Triangle_3_Ray_3_do_intersect.h.

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_TRIANGLE_3_RAY_3_INTERSECTION_H
#define CGAL_INTERNAL_INTERSECTIONS_3_TRIANGLE_3_RAY_3_INTERSECTION_H

#include <CGAL/kernel_basic.h>
#include <CGAL/intersections.h>

namespace CGAL {

namespace Intersections {

namespace internal {

template <class K>
typename K::Point_3
t3r3_intersection_coplanar_aux(const typename K::Point_3& p,
                               const typename K::Vector_3& v,
                               const typename K::Point_3& a,
                               const typename K::Point_3& b,
                               const K& k)
{
  // Returns the intersection point between line (p,v) and line (a,b)
  //
  // preconditions:
  //   + p,v,a,b are coplanar

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

  const Vector_3 ab = vector(a,b);
  const Vector_3 pa = vector(p,a);

  const Vector_3 pa_ab = cross_product(pa,ab);
  const Vector_3 v_ab = cross_product(v,ab);

  const FT t = scalar_product(pa_ab,v_ab) / sq_length(v_ab);

  return ( p + t*v );
}


template <class K>
typename Intersection_traits<K, typename K::Triangle_3, typename K::Ray_3>::result_type
t3r3_intersection_coplanar_aux(const typename K::Point_3& a,
                               const typename K::Point_3& b,
                               const typename K::Point_3& c,
                               const typename K::Ray_3& r,
                               const bool negative_side,
                               const K& k)
{
  // This function is designed to clip r into the triangle abc.
  // Point configuration should be as follows
  //
  //
  //    p+    +b
  //     |
  //  +c |       +a
  //     |
  //     |r
  //
  // We know that c is isolated on the negative side of r
  // but we don't know p position wrt [bc]&[ca]

  typedef typename K::Point_3 Point_3;
  typedef typename K::Vector_3 Vector_3;

  typename K::Coplanar_orientation_3 coplanar_orientation =
    k.coplanar_orientation_3_object();

  typename K::Construct_segment_3 segment =
    k.construct_segment_3_object();

  typename K::Construct_point_on_3 point_on =
    k.construct_point_on_3_object();

  const Point_3& p = point_on(r,0);

  // A ray is not symetric, 2 cases depending on isolated side of c
  Orientation cap = negative_side ? coplanar_orientation(c,a,p)
                                  : coplanar_orientation(b,c,p);

  switch ( cap ) {

    case NEGATIVE:
      // p is bellow [c,a]
      return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>();

    case COLLINEAR:
      return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>(p);

    case POSITIVE:
    {
      // Compute the intersection points between ray and [b,c],[c,a]
      Vector_3 v = r.to_vector();

      // Get intersection point at p side
      Point_3 p_side_end_point(p);
      Point_3 q_side_end_point;

      // A ray is not symetric, 2 cases depending on isolated side of c
      if ( negative_side )
      {
        if ( NEGATIVE == coplanar_orientation(b,c,p) )
          p_side_end_point = t3r3_intersection_coplanar_aux(p,v,b,c,k);

        // Get other end point (always intersection computation on the unbounded
        // side of the ray)
        q_side_end_point = t3r3_intersection_coplanar_aux(p,v,c,a,k);
      }
      else
      {
        if ( NEGATIVE == coplanar_orientation(c,a,p) )
          p_side_end_point = t3r3_intersection_coplanar_aux(p,v,c,a,k);

        // Get other end point (always intersection computation on the unbounded
        // side of the ray)
        q_side_end_point = t3r3_intersection_coplanar_aux(p,v,b,c,k);
      }

      // Build result
      return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>(segment(p_side_end_point, q_side_end_point));
    }

    default: // should not happen.
      CGAL_error();
      return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>();
  }

  CGAL_error();
  return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>();
}

template <class K>
typename Intersection_traits<K, typename K::Triangle_3, typename K::Ray_3>::result_type
intersection_coplanar(const typename K::Triangle_3 &t,
                      const typename K::Ray_3  &r,
                      const K & k )
{
  CGAL_kernel_precondition( ! k.is_degenerate_3_object()(t) ) ;
  CGAL_kernel_precondition( ! k.is_degenerate_3_object()(r) ) ;

  typedef typename K::Point_3 Point_3;

  typename K::Construct_point_on_3 point_on =
    k.construct_point_on_3_object();

  typename K::Construct_vertex_3 vertex_on =
    k.construct_vertex_3_object();

  typename K::Coplanar_orientation_3 coplanar_orientation =
    k.coplanar_orientation_3_object();

  typename K::Construct_segment_3 segment =
    k.construct_segment_3_object();

  typename K::Collinear_are_ordered_along_line_3 collinear_ordered =
    k.collinear_are_ordered_along_line_3_object();


  const Point_3 & p = point_on(r,0);
  const Point_3 & q = point_on(r,1);

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

  // Test whether the ray's supporting line intersects the
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
              return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>();

            case NEGATIVE:
              // c is isolated on the negative side
              return t3r3_intersection_coplanar_aux(a,b,c,r,true,k);

            default: // COLLINEAR
              // p,q,c are collinear
              if ( collinear_ordered(p,c,q) || collinear_ordered(p,q,c) )
                return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>(c);
              else
                return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>();
          }

        case NEGATIVE:
          if ( POSITIVE == pqc )
            // b is isolated on the negative side
            return t3r3_intersection_coplanar_aux(c,a,b,r,true,k);
          else
            // a is isolated on the positive side (here mb c could be use as
            // an endpoint instead of computing an intersection is some cases)
            return t3r3_intersection_coplanar_aux(b,c,a,r,false,k);

        case COLLINEAR:
          switch ( pqc ) {
            case POSITIVE:
              // p,q,b are collinear
              if ( collinear_ordered(p,b,q) || collinear_ordered(p,q,b) )
                return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>(b);
              else
                return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>();

            case NEGATIVE:
              // a is isolated on the positive side (here mb b could be use as
              // an endpoint instead of computing an intersection)
              return t3r3_intersection_coplanar_aux(b,c,a,r,false,k);

            default: // COLLINEAR
              // b,c,p,q are aligned, [p,q]&[b,c] have the same direction
              if ( collinear_ordered(p,b,c) )
                return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>(segment(b,c));
              else
                return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>(segment(p,c));
          }

        default: // should not happen.
          CGAL_error();
          return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>();
      }

    // -----------------------------------
    // pqa NEGATIVE
    // -----------------------------------
    case NEGATIVE:
      switch ( pqb ) {
        case POSITIVE:
          if ( POSITIVE == pqc )
            // a is isolated on the negative side
            return t3r3_intersection_coplanar_aux(b,c,a,r,true,k);
          else
            // b is isolated on the positive side (here mb c could be use as
            // an endpoint instead of computing an intersection, in some cases)
            return t3r3_intersection_coplanar_aux(c,a,b,r,false,k);

        case NEGATIVE:
          switch ( pqc ) {
            case POSITIVE:
              // c is isolated on the positive side
              return t3r3_intersection_coplanar_aux(a,b,c,r,false,k);

            case NEGATIVE:
              // the triangle lies in the negative halfspace
              // defined by the segment's supporting line.
              return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>();

            default: // COLLINEAR:
              // p,q,c are collinear
              if ( collinear_ordered(p,c,q) || collinear_ordered(p,q,c) )
                return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>(c);
              else
                return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>();
          }

        case COLLINEAR:
          switch ( pqc ) {
            case POSITIVE:
              // a is isolated on the negative side (here mb b could be use as
              // an endpoint instead of computing an intersection)
              return t3r3_intersection_coplanar_aux(b,c,a,r,true,k);

            case NEGATIVE:
              // p,q,b are collinear
              if ( collinear_ordered(p,b,q) || collinear_ordered(p,q,b) )
                return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>(b);
              else
                return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>();

            default: //  COLLINEAR
              // b,c,p,q are aligned, [p,q]&[c,b] have the same direction
              if ( collinear_ordered(p,c,b) )
                return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>(segment(c,b));
              else
                return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>(segment(p,b));
          }

        default: // should not happen.
          CGAL_error();
          return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>();
      }

    // -----------------------------------
    // pqa COLLINEAR
    // -----------------------------------
    case COLLINEAR:
      switch ( pqb ) {
        case POSITIVE:
          switch ( pqc ) {
            case POSITIVE:
              // p,q,a are collinear
              if ( collinear_ordered(p,a,q) || collinear_ordered(p,q,a) )
                return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>(a);
              else
                return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>();

            case NEGATIVE:
              // b is isolated on the positive side (here mb a could be use as
              // an endpoint instead of computing an intersection)
              return t3r3_intersection_coplanar_aux(c,a,b,r,false,k);

            default: // COLLINEAR
              // a,c,p,q are aligned, [p,q]&[c,a] have the same direction
              if ( collinear_ordered(p,c,a) )
                return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>(segment(c,a));
              else
                return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>(segment(p,a));
          }

        case NEGATIVE:
          switch ( pqc ) {
            case POSITIVE:
              // b is isolated on the negative side (here mb a could be use as
              // an endpoint instead of computing an intersection)
              return t3r3_intersection_coplanar_aux(c,a,b,r,true,k);

            case NEGATIVE:
              // p,q,a are collinear
              if ( collinear_ordered(p,a,q) || collinear_ordered(p,q,a) )
                return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>(a);
              else
                return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>();

            default: // COLLINEAR
              // a,c,p,q are aligned, [p,q]&[a,c] have the same direction
              if ( collinear_ordered(p,a,c) )
                return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>(segment(a,c));
              else
                return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>(segment(p,c));
          }

        case COLLINEAR:
          switch ( pqc ) {
            case POSITIVE:
              // a,b,p,q are aligned, [p,q]&[a,b] have the same direction
              if ( collinear_ordered(p,a,b) )
                return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>(segment(a,b));
              else
                return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>(segment(p,b));

            case NEGATIVE:
              // a,b,p,q are aligned, [p,q]&[b,a] have the same direction
              if ( collinear_ordered(p,b,a) )
                return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>(segment(b,a));
              else
                return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>(segment(p,a));

            default: // COLLINEAR
              // case pqc == COLLINEAR is impossible since the triangle is
              // assumed to be non flat
              CGAL_error();
              return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>();
          }

        default: // should not happen.
          CGAL_error();
          return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>();
      }

    default:// should not happen.
      CGAL_error();
      return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>();
  }
}


template <class K>
inline
typename boost::optional<typename K::Point_3>
t3r3_intersection_aux(const typename K::Triangle_3 &t,
                      const typename K::Ray_3 &r,
                      const K& k)
{
  typename Intersection_traits<K, typename K::Line_3, typename K::Plane_3>
    ::result_type
    v = internal::intersection(r.supporting_line(),t.supporting_plane(), k);

  if(v) {
    if(const typename K::Point_3* p = intersect_get<typename K::Point_3>(v))
      return *p;
  }
  return boost::optional<typename K::Point_3>();
}


template <class K>
typename Intersection_traits<K, typename K::Triangle_3,
                               typename K::Ray_3>::result_type
intersection(const typename K::Triangle_3  &t,
             const typename K::Ray_3 &r,
             const K& k)
{
  CGAL_kernel_precondition( ! k.is_degenerate_3_object()(t) ) ;
  CGAL_kernel_precondition( ! k.is_degenerate_3_object()(r) ) ;

  typedef typename K::Point_3 Point_3;

  typename K::Construct_vertex_3 vertex_on =
    k.construct_vertex_3_object();

  typename K::Orientation_3 orientation =
    k.orientation_3_object();

  typename K::Construct_ray_3 ray =
    k.construct_ray_3_object();

  typename K::Construct_point_on_3 point_on =
    k.construct_point_on_3_object();


  const Point_3& a = vertex_on(t,0);
  const Point_3& b = vertex_on(t,1);
  const Point_3& c = vertex_on(t,2);
  const Point_3& p = point_on(r,0);
  const Point_3& q = point_on(r,1);

  Point_3 d = point_on(ray(a,r.to_vector()),1);

  const Orientation ray_direction = orientation(a,b,c,d);
  const Orientation abcp = orientation(a,b,c,p);

  switch ( abcp ) {
    case POSITIVE:
      switch ( ray_direction ) {
        case POSITIVE:
          // the ray lies in the positive open halfspaces defined by the
          // triangle's supporting plane
          return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>();

        case NEGATIVE:
          // The ray straddles the triangle's plane
          // p sees the triangle in counterclockwise order
          if ( orientation(p,q,a,b) != POSITIVE
              && orientation(p,q,b,c) != POSITIVE
               && orientation(p,q,c,a) != POSITIVE ) {
            boost::optional<Point_3> op = t3r3_intersection_aux(t,r,k);
            if(op) return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>(*op);
            else return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>();
          }
          else
            return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>();

        default: // COPLANAR
          // The ray lie in a plane parallel to a,b,c support plane
          return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>();
      }

    case NEGATIVE:
      switch ( ray_direction ) {
        case POSITIVE:
          // The ray straddles the triangle's plane
          // q sees the triangle in counterclockwise order

          if ( orientation(q,p,a,b) != POSITIVE
              && orientation(q,p,b,c) != POSITIVE
               && orientation(q,p,c,a) != POSITIVE ) {
            boost::optional<Point_3> op = t3r3_intersection_aux(t,r,k);
            if(op) return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>(*op);
            else return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>();
          }
          else
            return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>();

        case NEGATIVE:
          // the ray lies in the negative open halfspaces defined by the
          // triangle's supporting plane
          return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>();

        default: // COPLANAR
          // The ray lie in a plane parallel to a,b,c support plane
          return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>();

      }

    case COPLANAR: // p belongs to the triangle's supporting plane
      switch ( ray_direction ) {
        case POSITIVE:
          // q sees the triangle in counterclockwise order
          if ( orientation(q,p,a,b) != POSITIVE
              && orientation(q,p,b,c) != POSITIVE
              && orientation(q,p,c,a) != POSITIVE )
          {
            boost::optional<Point_3> op = t3r3_intersection_aux(t,r,k);
            if(op) return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>(*op);
            else return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>();
          }
          else
            return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>();

        case NEGATIVE:
          // q sees the triangle in clockwise order
          if ( orientation(p,q,a,b) != POSITIVE
              && orientation(p,q,b,c) != POSITIVE
              && orientation(p,q,c,a) != POSITIVE )
          {
            boost::optional<Point_3> op = t3r3_intersection_aux(t,r,k);
            if(op) return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>(*op);
            else return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>();
          }
          else
            return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>();

        default: //  COPLANAR
          // The ray lies in the triangle supporting plane
          return intersection_coplanar(t,r,k);
      }

    default: // should not happen.
      CGAL_error();
      return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>();
  }

  CGAL_error();
  return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Ray_3>();
}

template <class K>
typename Intersection_traits<K, typename K::Ray_3, typename K::Triangle_3>::result_type
intersection(const typename K::Ray_3 &r,
             const typename K::Triangle_3  &t,
             const K& k) {
  return Intersections::internal::intersection(t, r, k);
}

} // end namespace internal
} // namespace Intersections
} // end namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_TRIANGLE_3_RAY_3_INTERSECTION_H
