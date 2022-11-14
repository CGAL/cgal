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

#ifndef CGAL_INTERNAL_INTERSECTIONS_SEGMENT_3_TRIANGLE_3_DO_INTERSECT_H
#define CGAL_INTERNAL_INTERSECTIONS_SEGMENT_3_TRIANGLE_3_DO_INTERSECT_H

#include <CGAL/enum.h>
#include <CGAL/kernel_assertions.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
bool do_intersect_coplanar(const typename K::Point_3& A,
                           const typename K::Point_3& B,
                           const typename K::Point_3& C,
                           const typename K::Point_3& p,
                           const typename K::Point_3& q, const K& k)
{
  typedef typename K::Point_3 Point_3;

  typename K::Coplanar_orientation_3 coplanar_orientation = k.coplanar_orientation_3_object();

  const Point_3* a = &A;
  const Point_3* b = &B;
  const Point_3* c = &C;

  // Determine the orientation of the triangle in the common plane

  if (coplanar_orientation(A,B,C) != POSITIVE)
  {
    // The triangle is not counterclockwise oriented
    // swap two vertices.
    b = &C;
    c = &B;
  }

  // Test whether the segment's supporting line intersects the
  // triangle in the common plane

  const Orientation pqa = coplanar_orientation(p,q,*a);
  const Orientation pqb = coplanar_orientation(p,q,*b);
  const Orientation pqc = coplanar_orientation(p,q,*c);

  switch ( pqa ) {
  case POSITIVE:
    switch ( pqb ) {
    case POSITIVE:
      if (pqc == POSITIVE)
        return false;

        // the triangle lies in the positive halfspace
        // defined by the segment's supporting line.
      // c is isolated on the negative side
      return coplanar_orientation(*b,*c,q) != NEGATIVE
          && coplanar_orientation(*c,*a,p) != NEGATIVE;
    case NEGATIVE:
      if (pqc == POSITIVE) // b is isolated on the negative side
        return coplanar_orientation(*a,*b,q) != NEGATIVE
            && coplanar_orientation(*b,*c,p) != NEGATIVE;
      // a is isolated on the positive side
      return coplanar_orientation(*a,*b,q) != NEGATIVE
          && coplanar_orientation(*c,*a,p) != NEGATIVE;
    case COLLINEAR:
      if (pqc == POSITIVE) // b is isolated on the negative side
        return coplanar_orientation(*a,*b,q) != NEGATIVE
            && coplanar_orientation(*b,*c,p) != NEGATIVE;
      // a is isolated on the positive side
      return coplanar_orientation(*a,*b,q) != NEGATIVE
          && coplanar_orientation(*c,*a,p) != NEGATIVE;
    default:// should not happen.
      CGAL_kernel_assertion(false);
      return false;
    }
  case NEGATIVE:
    switch ( pqb ) {
    case POSITIVE:
      if (pqc == POSITIVE) // a is isolated on the negative side
        return coplanar_orientation(*a,*b,p) != NEGATIVE
            && coplanar_orientation(*c,*a,q) != NEGATIVE;
      // b is isolated on the positive side
      return coplanar_orientation(*a,*b,p) != NEGATIVE
          && coplanar_orientation(*b,*c,q) != NEGATIVE;
    case NEGATIVE:
      if (pqc == NEGATIVE)
        return false;

      // the triangle lies in the negative halfspace
      // defined by the segment's supporting line.
      // c is isolated on the positive side
      return coplanar_orientation(*b,*c,p) != NEGATIVE
          && coplanar_orientation(*c,*a,q) != NEGATIVE;
    case COLLINEAR:
      if (pqc == NEGATIVE) // b is isolated on the positive side
        return coplanar_orientation(*a,*b,p) != NEGATIVE
            && coplanar_orientation(*b,*c,q) != NEGATIVE;
      // a is isolated on the negative side
      return coplanar_orientation(*a,*b,p) != NEGATIVE
          && coplanar_orientation(*c,*a,q) != NEGATIVE;

    default:// should not happen.
      CGAL_kernel_assertion(false);
      return false;
    }
  case COLLINEAR:
    switch ( pqb ) {
    case POSITIVE:
      if (pqc == POSITIVE) // a is isolated on the negative side
        return coplanar_orientation(*a,*b,p) != NEGATIVE
            && coplanar_orientation(*c,*a,q) != NEGATIVE;
      // b is isolated on the positive side
      return coplanar_orientation(*a,*b,p) != NEGATIVE
          && coplanar_orientation(*b,*c,q) != NEGATIVE;
    case NEGATIVE:
      if (pqc == NEGATIVE) // a is isolated on the positive side
        return coplanar_orientation(*a,*b,q) != NEGATIVE
            && coplanar_orientation(*c,*a,p) != NEGATIVE;
      // b is isolated on the negative side
      return coplanar_orientation(*a,*b,q) != NEGATIVE
          && coplanar_orientation(*b,*c,p) != NEGATIVE;
    case COLLINEAR:
      if (pqc == POSITIVE) // c is isolated on the positive side
        return coplanar_orientation(*b,*c,p) != NEGATIVE
            && coplanar_orientation(*c,*a,q) != NEGATIVE;
      // c is isolated on the negative side
      return coplanar_orientation(*b,*c,q) != NEGATIVE
           && coplanar_orientation(*c,*a,p) != NEGATIVE;
      // case pqc == COLLINEAR is impossible since the triangle is
      // assumed to be non flat

    default:// should not happen.
      CGAL_kernel_assertion(false);
      return false;

    }
  default:// should not happen.
    CGAL_kernel_assertion(false);
    return false;
  }
}

template <class K>
bool do_intersect_coplanar(const typename K::Triangle_3& t,
                           const typename K::Segment_3  &s,
                           const K& k )
{
  CGAL_kernel_precondition(!k.is_degenerate_3_object()(t));
  CGAL_kernel_precondition(!k.is_degenerate_3_object()(s));

  typedef typename K::Point_3 Point_3;

  typename K::Construct_point_on_3 point_on = k.construct_point_on_3_object();
  typename K::Construct_vertex_3 vertex_on = k.construct_vertex_3_object();

  const Point_3& p = point_on(s,0);
  const Point_3& q = point_on(s,1);

  const Point_3& A = vertex_on(t,0);
  const Point_3& B = vertex_on(t,1);
  const Point_3& C = vertex_on(t,2);

  return do_intersect_coplanar(A,B,C,p,q,k);
}

template <class K>
bool do_intersect(const typename K::Triangle_3& t,
                  const typename K::Segment_3& s,
                  const K& k)
{
  CGAL_kernel_precondition(!k.is_degenerate_3_object()(t) );
  CGAL_kernel_precondition(!k.is_degenerate_3_object()(s) );

  typedef typename K::Point_3 Point_3;

  typename K::Construct_point_on_3 point_on = k.construct_point_on_3_object();
  typename K::Construct_vertex_3 vertex_on = k.construct_vertex_3_object();
  typename K::Orientation_3 orientation = k.orientation_3_object();

  const Point_3& a = vertex_on(t,0);
  const Point_3& b = vertex_on(t,1);
  const Point_3& c = vertex_on(t,2);
  const Point_3& p = point_on(s,0);
  const Point_3& q = point_on(s,1);

  const Orientation abcp = orientation(a,b,c,p);
  const Orientation abcq = orientation(a,b,c,q);

  switch ( abcp ) {
  case POSITIVE:
    switch ( abcq ) {
    case POSITIVE:
      // the segment lies in the positive open halfspaces defined by the
      // triangle's supporting plane
      return false;
    case NEGATIVE:
      // p sees the triangle in counterclockwise order
      return orientation(p,q,a,b) != POSITIVE
          && orientation(p,q,b,c) != POSITIVE
          && orientation(p,q,c,a) != POSITIVE;
    case COPLANAR:
      // q belongs to the triangle's supporting plane
      // p sees the triangle in counterclockwise order
      return orientation(p,q,a,b) != POSITIVE
          && orientation(p,q,b,c) != POSITIVE
          && orientation(p,q,c,a) != POSITIVE;
    default: // should not happen.
      CGAL_kernel_assertion(false);
      return false;
    }
  case NEGATIVE:
    switch ( abcq ) {
    case POSITIVE:
      // q sees the triangle in counterclockwise order
      return orientation(q,p,a,b) != POSITIVE
          && orientation(q,p,b,c) != POSITIVE
          && orientation(q,p,c,a) != POSITIVE;
    case NEGATIVE:
      // the segment lies in the negative open halfspaces defined by the
      // triangle's supporting plane
      return false;
    case COPLANAR:
      // q belongs to the triangle's supporting plane
      // p sees the triangle in clockwise order
      return orientation(q,p,a,b) != POSITIVE
          && orientation(q,p,b,c) != POSITIVE
          && orientation(q,p,c,a) != POSITIVE;

    default: // should not happen.
      CGAL_kernel_assertion(false);
      return false;
    }
  case COPLANAR: // p belongs to the triangle's supporting plane
    switch ( abcq ) {
    case POSITIVE:
      // q sees the triangle in counterclockwise order
      return orientation(q,p,a,b) != POSITIVE
          && orientation(q,p,b,c) != POSITIVE
          && orientation(q,p,c,a) != POSITIVE;
    case NEGATIVE:
      // q sees the triangle in clockwise order
      return orientation(p,q,a,b) != POSITIVE
          && orientation(p,q,b,c) != POSITIVE
          && orientation(p,q,c,a) != POSITIVE;
    case COPLANAR:
      // the segment is coplanar with the triangle's supporting plane
      // we test whether the segment intersects the triangle in the common
      // supporting plane
      return do_intersect_coplanar(t,s,k);

    default: // should not happen.
      CGAL_kernel_assertion(false);
      return false;
    }
  default: // should not happen.
    CGAL_kernel_assertion(false);
    return false;
  }
}

template <class K>
inline
bool do_intersect(const typename K::Segment_3& s,
                  const typename K::Triangle_3& t,
                  const K& k)
{
  return do_intersect(t, s, k);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif //CGAL_INTERNAL_INTERSECTIONS_SEGMENT_3_TRIANGLE_3_DO_INTERSECT_H
