// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
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
//
//
// Author(s)     : Camille Wormser, Stephane Tayeb, Pierre Alliez
//

#ifndef CGAL_KERNEL_NEAREST_POINT_SEGMENT_3_H_
#define CGAL_KERNEL_NEAREST_POINT_SEGMENT_3_H_

#include <CGAL/kernel_basic.h>
#include <CGAL/enum.h>


namespace CGAL {


namespace internal {

/**
* @brief returns true if p is inside segment s. If p is not inside s,
* result is the nearest point of s from p. WARNING: it is assumed that
* t and p are on the same line.
* @param query the query point
* @param s the segment
* @param closest_point_on_segment if query is not inside s, the nearest point of s from p
* @param k the kernel
* @return true if p is inside s
*/
template <class K>
inline
bool
is_inside_segment_3(const typename K::Point_3& query,
                    const typename K::Segment_3 & s,
                    typename K::Point_3& closest_point_on_segment,
                    const K& k)
{
  typename K::Construct_vector_3 vector =
    k.construct_vector_3_object();
  typename K::Construct_vertex_3 vertex_on =
    k.construct_vertex_3_object();
  typename K::Compute_scalar_product_3 scalar_product =
    k.compute_scalar_product_3_object();

  typedef typename K::FT FT;
  typedef typename K::Point_3 Point;

  const Point& a = vertex_on(s, 0);
  const Point& b = vertex_on(s, 1);
  if( scalar_product(vector(a,b), vector(a, query)) < FT(0) )
  {
    closest_point_on_segment = a;
    return false;
  }
  if( scalar_product(vector(b,a), vector(b, query)) < FT(0) )
  {
    closest_point_on_segment = b;
    return false;
  }

  // query is on segment
  return true;
}

template <class K>
typename K::Point_3
nearest_point_3(const typename K::Point_3& query,
                const typename K::Segment_3& segment,
                const K& k)
{
  typedef typename K::Point_3 Point_3;

  typename K::Construct_projected_point_3 projection =
      k.construct_projected_point_3_object();
  typename K::Is_degenerate_3 is_degenerate =
      k.is_degenerate_3_object();
  typename K::Construct_vertex_3 vertex =
      k.construct_vertex_3_object();

  if(is_degenerate(segment))
    return vertex(segment, 0);

  // Project query on segment supporting line
  const Point_3 proj = projection(segment.supporting_line(), query);

  Point_3 closest_point_on_segment;
  bool inside = is_inside_segment_3(proj,segment,closest_point_on_segment,k);

  // If proj is inside segment, returns it
  if ( inside )
    return proj;

  // Else returns the constructed point
  return closest_point_on_segment;
}

} }  // end namespace CGAL::internal

#endif // CGAL_KERNEL_NEAREST_POINT_SEGMENT_3_H_
