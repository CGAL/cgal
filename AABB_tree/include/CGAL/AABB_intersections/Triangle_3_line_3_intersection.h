// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
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
// $URL: $
// $Id: $
//
//
// Author(s)     : Stéphane Tayeb
//
//******************************************************************************
// File Description :
//
//******************************************************************************

#ifndef TRIANGLE_3_LINE_3_INTERSECTION_H_
#define TRIANGLE_3_LINE_3_INTERSECTION_H_

#include <CGAL/kernel_basic.h>
#include <CGAL/intersections.h>

namespace CGAL {
namespace CGALi {

template <class K>
Object
intersection(const typename K::Triangle_3 &t,
             const typename K::Line_3 &l,
             const K& k)
{
  typedef typename K::Point_3 Point_3;
  typedef typename K::Segment_3 Segment_3;
  typedef typename K::Object_3 Object_3;
  typedef typename K::Line_3 Line_3;

  if ( !CGAL::do_intersect(t, l) )
    return Object_3();

  // TOFIX: here we assume that we have already tested
  // do_intersect between the triangle and the line
  const Object intersection =  CGAL::intersection(t.supporting_plane(),
                                                  l);

  // intersection is either a point, either a line
  // if it is a line, then we need to clip it to a segment
  if ( const Line_3* line = object_cast<Line_3>(&intersection))
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

    const Object_3 inter_01 = intersect(l01,*line);
    const Object_3 inter_02 = intersect(l02,*line);
    const Object_3 inter_12 = intersect(l12,*line);

    const Point_3* p01 = object_cast<Point_3>(&inter_01);
    const Point_3* p02 = object_cast<Point_3>(&inter_02);
    const Point_3* p12 = object_cast<Point_3>(&inter_12);

    if ( p01 && has_on(s01, *p01) )
    {
      if ( p02 && has_on(s02, *p02) )
        return make_object(f_segment(*p01,*p02));
      else if ( p12 && has_on(s12, *p12) )
        return make_object(f_segment(*p01,*p12));
      else
        return make_object(*p01);
    }
    else if ( p02 && has_on(s02, *p02) )
    {
      if  ( p12 && has_on(s12, *p12) )
        return make_object(f_segment(*p02,*p12));
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
intersection(const Triangle_3<K> &t, const Line_3<K> &l)
{
  return typename K::Intersect_3()(t,l);
}

template <class K>
inline
Object
intersection(const Line_3<K> &l, const Triangle_3<K> &t)
{
  return typename K::Intersect_3()(t,l);
}

} // end namespace CGAL

#endif // TRIANGLE_3_LINE_3_INTERSECTION_H_
