// Copyrigth (c) 2009  GeometryFactory (France)
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

#ifndef CGAL_TRIANGLE_3_RAY_3_INTERSECTION_H
#define CGAL_TRIANGLE_3_RAY_3_INTERSECTION_H

#include <CGAL/kernel_basic.h>
#include <CGAL/intersections.h>

namespace CGAL {
namespace CGALi {
template <class K>
Object
intersection(const typename K::Triangle_3  &t,
             const typename K::Ray_3 &r,
             const K& k)
{
  typedef typename K::Point_3 Point_3;
  typedef typename K::Segment_3 Segment_3;
  typedef typename K::Object_3 Object_3;
  typedef typename K::Ray_3 Ray_3;
  typedef typename K::Line_3 Line_3;


  if ( !CGAL::do_intersect(t, r) )
    return Object_3();

	// TOFIX: here we assume that we have already tested
	// do_intersect between the triangle and the ray
  const Object intersection =  CGAL::intersection(t.supporting_plane(),
                                                  r);

  // intersection is either a point, either a ray
  // if it is a ray, then we need to clip it to a segment
  if ( const Ray_3* ray = object_cast<Ray_3>(&intersection))
  {
    typename K::Intersect_3 intersect = k.intersect_3_object();
    typename K::Construct_line_3 line = k.construct_line_3_object();
    typename K::Construct_vertex_3 vertex_on = k.construct_vertex_3_object();
    typename K::Construct_object_3 make_object = k.construct_object_3_object();
    typename K::Construct_segment_3 segment = k.construct_segment_3_object();

    const Point_3& t0 = vertex_on(t,0);
    const Point_3& t1 = vertex_on(t,1);
    const Point_3& t2 = vertex_on(t,2);

    const Line_3 l01 = line(t0,t1);
    const Line_3 l02 = line(t0,t2);
    const Line_3 l12 = line(t1,t2);

    const Line_3 lr = ray->supporting_line();

    const Object_3 inter_01 = intersect(l01,lr);
    const Object_3 inter_02 = intersect(l02,lr);
    const Object_3 inter_12 = intersect(l12,lr);

    const Point_3* p01 = object_cast<Point_3>(&inter_01);
    const Point_3* p02 = object_cast<Point_3>(&inter_02);
    const Point_3* p12 = object_cast<Point_3>(&inter_12);

    if ( p01 )
    {
      if ( p02 )
        return make_object(segment(*p01,*p02));
      if ( p12 )
        return make_object(segment(*p01,*p12));
    }
    else
    {
      if ( p02 && p12)
        return make_object(segment(*p02,*p12));
    }

    // here ray intersects the boundary of the triangle
    if ( const Segment_3* s01 = object_cast<Segment_3>(&inter_01) )
      return make_object(*s01);
    if ( const Segment_3* s02 = object_cast<Segment_3>(&inter_02) )
      return make_object(*s02);
    if ( const Segment_3* s12 = object_cast<Segment_3>(&inter_12) )
      return make_object(*s12);

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
intersection(const Triangle_3<K> &t, const Ray_3<K> &r)
{
  return typename K::Intersect_3()(t, r);
}

template <class K>
inline
Object
intersection(const Ray_3<K> &r, const Triangle_3<K> &t)
{
  return typename K::Intersect_3()(t, r);
}

} // end namespace CGAL

#endif // CGAL_TRIANGLE_3_RAY_3_INTERSECTION_H
