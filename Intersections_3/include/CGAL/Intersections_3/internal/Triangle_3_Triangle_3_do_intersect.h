// Copyright (c) 2003  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Philippe Guigue

#ifndef CGAL_TRIANGLE_3_TRIANGLE_3_DO_INTERSECT_H
#define CGAL_TRIANGLE_3_TRIANGLE_3_DO_INTERSECT_H

#include <CGAL/Uncertain.h>
#include <CGAL/enum.h>
#include <CGAL/kernel_assertions.h>
#include <CGAL/Intersection_traits_3.h>

namespace CGAL {

  template <class K>
  class Triangle_3;

namespace Intersections {

namespace internal {

template <class K>
bool  _intersection_test_vertex(const typename K::Point_3 * p,
                                const typename K::Point_3 * q,
                                const typename K::Point_3 * r,
                                const typename K::Point_3 * a,
                                const typename K::Point_3 * b,
                                const typename K::Point_3 * c,
                                const K & k){

  CGAL_kernel_precondition( k.coplanar_orientation_3_object() (*p,*q,*r)
                            == POSITIVE);
  CGAL_kernel_precondition( k.coplanar_orientation_3_object() (*a,*b,*c)
                            == POSITIVE);


  // Vertex p sees vertex c

  typename K::Coplanar_orientation_3 coplanar_orientation =
    k.coplanar_orientation_3_object();


  if (coplanar_orientation(*c,*a,*q) != NEGATIVE) {
    if (coplanar_orientation(*c,*b,*q) != POSITIVE) {
      if (coplanar_orientation(*p,*a,*q) == POSITIVE)
        return  coplanar_orientation(*p,*b,*q) != POSITIVE;

      return coplanar_orientation(*p,*a,*r) != NEGATIVE
        && coplanar_orientation(*q,*r,*a) != NEGATIVE;
    }

    if (coplanar_orientation(*p,*b,*q) != POSITIVE)
      return coplanar_orientation(*c,*b,*r) != POSITIVE
        && coplanar_orientation(*q,*r,*b) != NEGATIVE;
    return false;

  }

  if (coplanar_orientation(*c,*a,*r) != NEGATIVE) { //qr straddles (ac)
    if (coplanar_orientation(*q,*r,*c) != NEGATIVE)
      return (coplanar_orientation(*p,*a,*r) != NEGATIVE);

    return coplanar_orientation(*q,*r,*b) != NEGATIVE
      && coplanar_orientation(*c,*r,*b) != NEGATIVE;
  }
  return false; // ca separes

}


template <class K>
bool  _intersection_test_edge(const typename K::Point_3 * p,
                              const typename K::Point_3 * q,
                              const typename K::Point_3 * r,
                              const typename K::Point_3 * a,
                              const typename K::Point_3 * CGAL_kernel_precondition_code(b),
                              const typename K::Point_3 * c,
                              const K & k){

  CGAL_kernel_precondition( k.coplanar_orientation_3_object() (*p,*q,*r)
                            == POSITIVE);
  CGAL_kernel_precondition( k.coplanar_orientation_3_object() (*a,*b,*c)
                            == POSITIVE);

  // Vertex p sees edge ca

  typename K::Coplanar_orientation_3 coplanar_orientation =
    k.coplanar_orientation_3_object();



  if (coplanar_orientation(*c,*a,*q) != NEGATIVE) {  //pq straddles (ac)
    if (coplanar_orientation(*p,*a,*q) != NEGATIVE)
      return coplanar_orientation(*p,*q,*c) != NEGATIVE ;

    return coplanar_orientation(*q,*r,*a) != NEGATIVE
      && coplanar_orientation(*r,*p,*a) != NEGATIVE;
  }

  if (coplanar_orientation(*c,*a,*r) != NEGATIVE) {
    // pr and qr straddle line (pr)
    return coplanar_orientation(*p,*a,*r) != NEGATIVE
      && ( coplanar_orientation(*p,*r,*c) != NEGATIVE
           || coplanar_orientation(*q,*r,*c) != NEGATIVE );
  }

  return false; //ac separes
}


template <class K>
bool do_intersect_coplanar(const typename K::Triangle_3 &t1,
                           const typename K::Triangle_3 &t2,
                           const K & k)
{

  CGAL_kernel_precondition( ! k.is_degenerate_3_object() (t1) );
  CGAL_kernel_precondition( ! k.is_degenerate_3_object() (t2) );

  typedef typename K::Point_3 Point_3;

  typename K::Construct_vertex_3 vertex_on =
    k.construct_vertex_3_object();

  typename K::Coplanar_orientation_3 coplanar_orientation =
    k.coplanar_orientation_3_object();

    const Point_3 & P = vertex_on(t1,0);
    const Point_3 & Q = vertex_on(t1,1);
    const Point_3 & R = vertex_on(t1,2);

    const Point_3 & A = vertex_on(t2,0);
    const Point_3 & B = vertex_on(t2,1);
    const Point_3 & C = vertex_on(t2,2);


    const Point_3 * p = &P;
    const Point_3 * q = &Q;
    const Point_3 * r = &R;

    const Point_3 * a = &A;
    const Point_3 * b = &B;
    const Point_3 * c = &C;

    // First we ensure that both triangles are counterclockwise
    // oriented, t1 and t2 are supposed  to be non flat.

    if ( coplanar_orientation(P,Q,R) == NEGATIVE ) {
      q = &R;
      r = &Q;
    }

    if ( coplanar_orientation(A,B,C) == NEGATIVE ) {
      b = &C;
      c = &B;
    }


    // Localization of p in the arrangement of the
    // lines supporting abc's edges

    if ( coplanar_orientation(*a,*b,*p) != NEGATIVE ) {
      if ( coplanar_orientation(*b,*c,*p) != NEGATIVE ) {
        if ( coplanar_orientation(*c,*a,*p) != NEGATIVE )
          // p is inside triangle abc
          return true;
        // p sees ac
        return internal::_intersection_test_edge(p,q,r,a,b,c,k);
      }
      if ( coplanar_orientation(*c,*a,*p) != NEGATIVE )//p sees bc
        return internal::_intersection_test_edge(p,q,r,c,a,b,k);
      // p sees c
      return internal::_intersection_test_vertex(p,q,r,a,b,c,k);

    }
    if ( coplanar_orientation(*b,*c,*p) != NEGATIVE ) {
      if ( coplanar_orientation(*c,*a,*p) != NEGATIVE ) //p sees ab
        return internal::_intersection_test_edge(p,q,r,b,c,a,k);
      // p sees a
      return internal::_intersection_test_vertex(p,q,r,b,c,a,k);
    }
    // p sees b
    return internal::_intersection_test_vertex(p,q,r,c,a,b,k);

}


template <class K>
typename K::Boolean
do_intersect(const typename K::Triangle_3 &t1,
             const typename K::Triangle_3 &t2,
             const K & k)
{
  CGAL_kernel_precondition( ! k.is_degenerate_3_object() (t1) );
  CGAL_kernel_precondition( ! k.is_degenerate_3_object() (t2) );

  typedef typename K::Point_3 Point_3;

  typename K::Construct_vertex_3 vertex_on =
    k.construct_vertex_3_object();

  typename K::Orientation_3 orientation =
    k.orientation_3_object();


   const Point_3 & p = vertex_on(t1,0);
   const Point_3 & q = vertex_on(t1,1);
   const Point_3 & r = vertex_on(t1,2);
   const Point_3 & a = vertex_on(t2,0);
   const Point_3 & b = vertex_on(t2,1);
   const Point_3 & c = vertex_on(t2,2);


   const Point_3  * s_min1;
   const Point_3  * t_min1;
   const Point_3  * s_max1;
   const Point_3  * t_max1;

   // Compute distance signs  of p, q and r to the plane of triangle(a,b,c)
   const Orientation dp = make_certain(orientation(a,b,c,p));
   const Orientation dq = make_certain(orientation(a,b,c,q));
   const Orientation dr = make_certain(orientation(a,b,c,r));


    switch ( dp ) {
    case POSITIVE:
      if ( dq == POSITIVE ) {
        if ( dr == POSITIVE)
          // pqr lies in the open positive halfspace induced by
          // the plane of triangle(a,b,c)
          return false;
        // r is isolated on the negative side of the plane
        s_min1 = &q ; t_min1 = &r ; s_max1 = &r ; t_max1 = &p;

      } else {
        if  ( dr == POSITIVE) {
          // q is isolated on the negative side of the plane
          s_min1 =  &p; t_min1 =  &q; s_max1 =  &q; t_max1 =  &r;
        } else {
          // p is isolated on the positive side of the plane
          s_min1 =  &p; t_min1 =  &q; s_max1 =  &r; t_max1 =  &p;
        }
      }
      break;

    case NEGATIVE:
      if ( dq == NEGATIVE ) {
        if ( dr == NEGATIVE )
          // pqr lies in the open negative halfspace induced by
          // the plane of triangle(a,b,c)
          return false;
        // r is isolated on the positive side of the plane
        s_min1 =  &r; t_min1 =  &p; s_max1 =  &q; t_max1 =  &r;

      } else {
        if ( dr == NEGATIVE ) {
          // q is isolated on the positive side of the plane
          s_min1 =  &q; t_min1 =  &r; s_max1 =  &p; t_max1 =  &q;
        } else {
          // p is isolated on the negative side of the plane
          s_min1 =  &r; t_min1 =  &p; s_max1 =  &p; t_max1 =  &q;
        }
      }
      break;

    case COPLANAR:
      switch  ( dq ) {
      case POSITIVE:
        if ( dr == POSITIVE ) {
          // p is isolated on the negative side of the plane
          s_min1 =  &r; t_min1 =  &p; s_max1 =  &p; t_max1 =  &q;
        } else {
          // q is isolated on the positive side of the plane
          s_min1 =  &q; t_min1 =  &r; s_max1 =  &p; t_max1 =  &q;
        }
        break;

      case NEGATIVE:
        if ( dr == NEGATIVE ) {
          // p is isolated on the positive side of the plane
          s_min1 =  &p; t_min1 =  &q; s_max1 =  &r; t_max1 =  &p;
        } else {
          // q is isolated on the negative side of the plane
          s_min1 =  &p; t_min1 =  &q; s_max1 =  &q; t_max1 =  &r;
        }
        break;

      case COPLANAR:
        switch ( dr ) {
        case POSITIVE:
          // r is isolated on the positive side of the plane
          s_min1 =  &r; t_min1 =  &p; s_max1 =  &q; t_max1 =  &r;
          break;

        case NEGATIVE:
          // r is isolated on the negative side of the plane
          s_min1 =  &q ; t_min1 =  &r ; s_max1 =  &r ; t_max1 =  &p;
          break;

        case COPLANAR:
          return do_intersect_coplanar(t1,t2,k);
        default: // should not happen.
          CGAL_kernel_assertion(false);
          return false;

        }
        break;

      default: // should not happen.
      CGAL_kernel_assertion(false);
      return false;

      }
      break;

    default: // should not happen.
      CGAL_kernel_assertion(false);
      return false;
    }



    const Point_3  * s_min2;
    const Point_3  * t_min2;
    const Point_3  * s_max2;
    const Point_3  * t_max2;

    // Compute distance signs  of a, b and c to the plane of triangle(p,q,r)
    const Orientation da = make_certain(orientation(p,q,r,a));
    const Orientation db = make_certain(orientation(p,q,r,b));
    const Orientation dc = make_certain(orientation(p,q,r,c));


    switch ( da ) {
    case POSITIVE:
      if ( db == POSITIVE ) {
        if ( dc == POSITIVE)
          // abc lies in the open positive halfspace induced by
          // the plane of triangle(p,q,r)
          return false;
        // c is isolated on the negative side of the plane
        s_min2 =  &b ; t_min2 =  &c ; s_max2 =  &c ; t_max2 =  &a;

      } else {
        if  ( dc == POSITIVE) {
          // b is isolated on the negative side of the plane
          s_min2 =  &a; t_min2 =  &b; s_max2 =  &b; t_max2 =  &c;
        } else {
          // a is isolated on the positive side of the plane
          s_min2 =  &a; t_min2 =  &b; s_max2 =  &c; t_max2 =  &a;
        }
      }
      break;

    case NEGATIVE:
      if ( db == NEGATIVE ) {
        if ( dc == NEGATIVE )
          // abc lies in the open negative halfspace induced by
          // the plane of triangle(p,q,r)
          return false;
        // c is isolated on the positive side of the plane
        s_min2 =  &c ; t_min2 =  &a ; s_max2 =  &b ; t_max2 =  &c;

      } else {
        if ( dc == NEGATIVE ) {
          // b is isolated on the positive side of the plane
          s_min2 =  &b ; t_min2 =  &c ; s_max2 =  &a ; t_max2 =  &b;
        } else {
          // a is isolated on the negative side of the plane
          s_min2 =  &c ; t_min2 =  &a ; s_max2 =  &a ; t_max2 =  &b;
        }
      }
      break;

    case COPLANAR:
      switch  ( db ) {
      case POSITIVE:
        if ( dc == POSITIVE ) {
          // a is isolated on the negative side of the plane
          s_min2 =  &c ; t_min2 =  &a ; s_max2 =  &a ; t_max2 =  &b;
        } else {
          // b is isolated on the positive side of the plane
          s_min2 =  &b ; t_min2 =  &c ; s_max2 =  &a ; t_max2 =  &b;
        }
        break;

      case NEGATIVE:
        if ( dc == NEGATIVE ) {
          // a is isolated on the positive side of the plane
          s_min2 =  &a ; t_min2 =  &b ; s_max2 =  &c ; t_max2 =  &a;
        } else {
          // b is isolated on the negative side of the plane
          s_min2 =  &a ; t_min2 =  &b ; s_max2 =  &b ; t_max2 =  &c;
        }
        break;

      case COPLANAR:
        switch ( dc ) {
        case POSITIVE:
          // c is isolated on the positive side of the plane
          s_min2 =  &c ; t_min2 =  &a ; s_max2 =  &b ; t_max2 =  &c;
          break;

        case NEGATIVE:
          // c is isolated on the negative side of the plane
          s_min2 =  &b ; t_min2 =  &c ; s_max2 =  &c ; t_max2 =  &a;
          break;

        case COPLANAR:
          // Supposed unreachable code
          // since the triangles are assumed to be non-flat

          return do_intersect_coplanar(t1,t2,k);
        default: // should not happen.
          CGAL_kernel_assertion(false);
          return false;

        }
        break;

      default: // should not happen.
      CGAL_kernel_assertion(false);
      return false;

      }
      break;

    default: // should not happen.
      CGAL_kernel_assertion(false);
      return false;
    }

    return  orientation(*s_min1,*t_min1,*s_min2,*t_min2) != POSITIVE
      &&  orientation(*s_max1,*t_max1,*t_max2,*s_max2) != POSITIVE;
}

} // namespace internal
} // namespace Intersections
} //namespace CGAL

#endif // CGAL_TRIANGLE_3_TRIANGLE_3_DO_INTERSECT_H
