// Copyright (c) 2010,2012 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Laurent Rineau
//

#ifndef CGAL_SQUARED_DISTANCE_POINT_3_TRIANGLE_3_H
#define CGAL_SQUARED_DISTANCE_POINT_3_TRIANGLE_3_H

#include <CGAL/squared_distance_3_0.h>
#include <CGAL/squared_distance_3_1.h>
#include <CGAL/wmult.h>

#include <CGAL/Point_3.h>
#include <CGAL/Triangle_3.h>

namespace CGAL {
namespace internal {

template <class K>
inline bool
on_left_of_triangle_edge(const typename K::Point_3 & pt,
                         const typename K::Vector_3 & normal,
                         const typename K::Point_3 & ep0,
                         const typename K::Point_3 & ep1,
                         const K& k)
{
  // return true iff pt is on the negative side of the plane defined
  // by (ep0, ep1) and normal
  typename K::Construct_vector_3 vector;
  typename K::Vector_3 edge = vector(ep0, ep1);
  typename K::Vector_3 diff = vector(ep0, pt);

  typedef typename K::RT RT;

  const bool result = 
    RT(wdot(wcross(edge,
                   normal,
                   k), 
            diff,
            k)) <= RT(0);
  return result;
}

template <class K>
inline typename K::FT
squared_distance_to_triangle(
    const typename K::Point_3 & pt,
    const typename K::Point_3 & t0,
    const typename K::Point_3 & t1,
    const typename K::Point_3 & t2,
    const K& k)
{
  typename K::Construct_vector_3 vector;
  typedef typename K::Vector_3 Vector_3;
  const Vector_3 e1 = vector(t0, t1);
  const Vector_3 oe3 = vector(t0, t2);
  const Vector_3 normal = wcross(e1, oe3, k);

  if(normal != NULL_VECTOR
     && on_left_of_triangle_edge(pt, normal, t0, t1, k)
     && on_left_of_triangle_edge(pt, normal, t1, t2, k)
     && on_left_of_triangle_edge(pt, normal, t2, t0, k))
      {
        // the projection of pt is inside the triangle
        return squared_distance_to_plane(normal, vector(t0, pt), k);
      }
      else {
        // The case normal==NULL_VECTOR covers the case when the triangle
        // is colinear, or even more degenerate. In that case, we can
        // simply take also the distance to the three segments.
        typename K::FT d1 = squared_distance(pt, 
                                             typename K::Segment_3(t2, t0),
                                             k);
        typename K::FT d2 = squared_distance(pt, 
                                             typename K::Segment_3(t1, t2),
                                             k);
        typename K::FT d3 = squared_distance(pt, 
                                             typename K::Segment_3(t0, t1),
                                             k);
       
        return (std::min)( (std::min)(d1, d2), d3);
      }
}

template <class K>
inline typename K::FT
squared_distance(
    const typename K::Point_3 & pt,
    const typename K::Triangle_3 & t,
    const K& k)
{
  typename K::Construct_vertex_3 vertex;
  return squared_distance_to_triangle(pt,
                                      vertex(t, 0),
                                      vertex(t, 1),
                                      vertex(t, 2),
                                      k);
}

} // end namespace CGAL::internal 

template <class K>
inline typename K::FT
squared_distance(const Point_3<K> & pt,
                 const Triangle_3<K> & t) {
  return internal::squared_distance(pt, t, K());
}


template <class K>
inline typename K::FT
squared_distance(const Triangle_3<K> & t,
                 const Point_3<K> & pt) {
  return internal::squared_distance(pt, t, K());
}

} // end namespace CGAL

#endif // CGAL_SQUARED_DISTANCE_POINT_3_TRIANGLE_3_H
