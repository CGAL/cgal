// Copyright (c) 2009  GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     :  Laurent Rineau

#ifndef CGAL_TRIANGLE_3_SEGMENT_3_INTERSECTION_H
#define CGAL_TRIANGLE_3_SEGMENT_3_INTERSECTION_H

#include <CGAL/kernel_basic.h>
#include <CGAL/intersection_3.h>

namespace CGAL {
namespace CGALi {



template <class K>
inline
Object
intersection_triangle_segment_aux(const typename K::Segment_3 &s1,
                                  const typename K::Segment_3 &s2,
                                  const K& k)
{
  typedef typename K::Point_3 Point_3;
  typedef typename K::Object_3 Object_3;
  typedef typename K::FT FT;

  typename K::Construct_object_3 make_object = k.construct_object_3_object();
  typename K::Construct_segment_3 make_segment = k.construct_segment_3_object();
  typename K::Compute_squared_distance_3 sq_distance = k.compute_squared_distance_3_object();

  const Point_3& p1 = s1.source();
  const Point_3& p2 = s1.target();
  const Point_3& q1 = s2.source();
  const Point_3& q2 = s2.target();

  // distance test are better than coordinates comparisons
  // (I had problems with has_on and are_colinear_ordered_along_line_3,
  // because some points are constructed by intersection)
  const FT p1p2 = sq_distance(p1,p2);
  const FT p1q1 = sq_distance(p1,q1);
  const FT p2q1 = sq_distance(p2,q1);
  const FT p1q2 = sq_distance(p1,q2);
  const FT p2q2 = sq_distance(p2,q2);

  if ( p1p2<p1q1 || p1p2<p2q1 ) // q1 is outside [p1,p2]
  {
    if ( p1p2<p1q2 || p1p2<p2q2 ) // q2 is outside [p1,p2]
    {
      if ( p1q1<p2q1 && p1q2<p2q2 ) // q1&q2 are on p1 side
        return make_object(p1);
      else if ( p2q1<p1q1 && p2q2<p1q2 ) // q1&q2 are on p2 side
        return make_object(p2);
      else // [p1,p2] is inside [q1,q2]
        return make_object(make_segment(p1,p2));
    }
    else // q2 is inside [p1,p2]
    {
      if ( p1q1 < p2q1 ) // q1 is on p1 side
        return make_object(make_segment(p1,q2));
      else // q1 is on p2 side
        return make_object(make_segment(p2,q2));
    }
  }
  else // q1 is inside [p1,p2]
  {
    if ( p1p2<p1q2 || p1p2<p2q2) // q2 is outside [p1,p2]
    {
      if ( p1q2<p2q2 ) // q2 is on p1 side
        return make_object(make_segment(p1,q1));
      else // q2 is on p2 side
        return make_object(make_segment(p2,q1));
    }
    else // q2 is inside [p1,p2]
      return make_object(make_segment(q1,q2));
  }

  return Object(); // unreachable code
}


template <class K>
Object
intersection(const typename K::Triangle_3  &t,
	     const typename K::Segment_3 &s,
	     const K& k)
{
  typedef typename K::Point_3 Point_3;
  typedef typename K::Segment_3 Segment_3;
  typedef typename K::Object_3 Object_3;
  typedef typename K::Ray_3 Ray_3;
  typedef typename K::Line_3 Line_3;


  if ( !CGAL::do_intersect(t, s) )
    return Object_3();

  // TOFIX: here we assume that we have already tested
  // do_intersect between the triangle and the segment
  const Object intersection =  CGAL::intersection(t.supporting_plane(),
                                                  s);

  // intersection is either a point, either a segment
  // if it is a segment, then we need to clip it to a segment located inside
  // the triangle
  if ( const Segment_3* seg = object_cast<Segment_3>(&intersection))
  {
    typename K::Intersect_3 intersect = k.intersect_3_object();
    typename K::Construct_line_3 f_line = k.construct_line_3_object();
    typename K::Construct_vertex_3 vertex_on = k.construct_vertex_3_object();
    typename K::Construct_object_3 make_object = k.construct_object_3_object();
    typename K::Construct_segment_3 f_segment = k.construct_segment_3_object();
    typename K::Has_on_3 has_on = k.has_on_3_object();

    const Point_3& t0 = vertex_on(t,0);
    const Point_3& t1 = vertex_on(t,1);
    const Point_3& t2 = vertex_on(t,2);

    const Line_3 l01 = f_line(t0,t1);
    const Line_3 l02 = f_line(t0,t2);
    const Line_3 l12 = f_line(t1,t2);

    const Segment_3 s01 = f_segment(t0,t1);
    const Segment_3 s02 = f_segment(t0,t2);
    const Segment_3 s12 = f_segment(t1,t2);

    const Line_3 ls = seg->supporting_line();

    const Object_3 inter_01 = intersect(l01,ls);
    const Object_3 inter_02 = intersect(l02,ls);
    const Object_3 inter_12 = intersect(l12,ls);

    const Point_3* p01 = object_cast<Point_3>(&inter_01);
    const Point_3* p02 = object_cast<Point_3>(&inter_02);
    const Point_3* p12 = object_cast<Point_3>(&inter_12);

    // has_on is probably not robust (see below)
    if ( p01 && has_on(s01, *p01) )
    {
      if ( p02 && has_on(s02, *p02) )
        return intersection_triangle_segment_aux(*seg, f_segment(*p01,*p02), k);
      else if ( p12 && has_on(s12, *p12) )
        return intersection_triangle_segment_aux(*seg, f_segment(*p01,*p12), k);
      else
        return make_object(*p01);
    }
    else if ( p02 && has_on(s02, *p02) )
    {
      if  ( p12 && has_on(s12, *p12) )
        return intersection_triangle_segment_aux(*seg, f_segment(*p02,*p12), k);
      else
        return make_object(*p02);
    }
    else if ( p12 && has_on(s12, *p12) )
    {
      return make_object(*p12);
    }

    // should not happen
    CGAL_kernel_assertion(false);
    return Object_3();
  }
  else
    return intersection;
}

} // end namespace CGALi

template <class K>
inline
Object
intersection(const Triangle_3<K> &t, const Segment_3<K> &s)
{
  return typename K::Intersect_3()(t, s);
}

template <class K>
inline
Object
intersection(const Segment_3<K> &s, const Triangle_3<K> &t)
{
  return typename K::Intersect_3()(t, s);
}

} // end namespace CGAL

#endif // CGAL_TRIANGLE_3_SEGMENT_3_INTERSECTION_H

