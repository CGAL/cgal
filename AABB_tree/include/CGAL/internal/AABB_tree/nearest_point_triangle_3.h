// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Camille Wormser, Stephane Tayeb
//
//******************************************************************************
// File Description :
//
//******************************************************************************

#ifndef NEAREST_POINT_TRIANGLE_3_H_
#define NEAREST_POINT_TRIANGLE_3_H_

#include <CGAL/kernel_basic.h>
#include <CGAL/enum.h>


namespace CGAL {


namespace internal {


template <class K>
inline
bool
is_inside_triangle_3_aux(const typename K::Vector_3& w,
                         const typename K::Point_3& p1,
                         const typename K::Point_3& p2,
                         const typename K::Point_3& q,
                         typename K::Point_3& result,
                         bool& outside,
                         const K& k)
{
  typedef typename K::Vector_3 Vector_3;
  typedef typename K::FT FT;

  typename K::Construct_vector_3 vector =
    k.construct_vector_3_object();
  typename K::Construct_projected_point_3 projection =
    k.construct_projected_point_3_object();
  typename K::Construct_line_3 line =
    k.construct_line_3_object();
  typename K::Compute_scalar_product_3 scalar_product =
    k.compute_scalar_product_3_object();
  typename K::Construct_cross_product_vector_3 cross_product =
    k.construct_cross_product_vector_3_object();

  const Vector_3 v = cross_product(vector(p1,p2), vector(p1,q));
  if ( scalar_product(v,w) < FT(0))
  {
    if (   scalar_product(vector(p1,q), vector(p1,p2)) >= FT(0)
        && scalar_product(vector(p2,q), vector(p2,p1)) >= FT(0) )
    {
      result = projection(line(p1, p2), q);
      return true;
    }
    outside = true;
  }

  return false;
}


/**
 * Returns the nearest point of p1,p2,p3 from origin
 * @param origin the origin point
 * @param p1 the first point
 * @param p2 the second point
 * @param p3 the third point
 * @param k the kernel
 * @return the nearest point from origin
 */
template <class K>
inline
typename K::Point_3
nearest_point_3(const typename K::Point_3& origin,
                const typename K::Point_3& p1,
                const typename K::Point_3& p2,
                const typename K::Point_3& p3,
                const K& k)
{
  typedef typename K::FT FT;

  typename K::Compute_squared_distance_3 sq_distance =
    k.compute_squared_distance_3_object();

  const FT dist_origin_p1 = sq_distance(origin,p1);
  const FT dist_origin_p2 = sq_distance(origin,p2);
  const FT dist_origin_p3 = sq_distance(origin,p3);

  if (   dist_origin_p2 >= dist_origin_p1
      && dist_origin_p3 >= dist_origin_p1 )
  {
    return p1;
  }
  if ( dist_origin_p3 >= dist_origin_p2 )
  {
    return p2;
  }

  return p3;
}

/**
 * @brief returns true if p is inside triangle t. If p is not inside t,
 * result is the nearest point of t from p. WARNING: it is assumed that
 * t and p are on the same plane.
 * @param p the reference point
 * @param t the triangle
 * @param result if p is not inside t, the nearest point of t from p
 * @param k the kernel
 * @return true if p is inside t
 */
template <class K>
inline
bool
is_inside_triangle_3(const typename K::Point_3& p,
                     const typename K::Triangle_3& t,
                     typename K::Point_3& result,
                     const K& k)
{
  typedef typename K::Point_3 Point_3;
  typedef typename K::Vector_3 Vector_3;

  typename K::Construct_vector_3 vector =
    k.construct_vector_3_object();
  typename K::Construct_vertex_3 vertex_on =
    k.construct_vertex_3_object();
  typename K::Construct_cross_product_vector_3 cross_product =
    k.construct_cross_product_vector_3_object();

  const Point_3& t0 = vertex_on(t,0);
  const Point_3& t1 = vertex_on(t,1);
  const Point_3& t2 = vertex_on(t,2);

  Vector_3 w = cross_product(vector(t0,t1), vector(t1,t2));

  bool outside = false;
  if (   is_inside_triangle_3_aux(w, t0, t1, p, result, outside, k)
      || is_inside_triangle_3_aux(w, t1, t2, p, result, outside, k)
      || is_inside_triangle_3_aux(w, t2, t0, p, result, outside, k) )
  {
    return false;
  }

  if ( outside )
  {
    result = nearest_point_3(p,t0,t1,t2,k);
    return false;
  }
  else
  {
    return true;
  }
}

/**
 * @brief Computes the closest_point from origin between bound and
 * any point of triangle.
 * @param origin the origin point
 * @param triangle the triangle
 * @param bound the farthest point
 * @param k the kernel
 * @return nearest point: bound or a point inside triangle
 */
template <class K>
typename K::Point_3
nearest_point_3(const typename K::Point_3& origin,
                const typename K::Triangle_3& triangle,
                const typename K::Point_3& bound,
                const K& k)
{
  typedef typename K::Point_3 Point_3;
  typedef typename K::FT FT;

  typename K::Compute_squared_distance_3 sq_distance =
    k.compute_squared_distance_3_object();
  typename K::Compare_squared_distance_3 compare_sq_distance =
    k.compare_squared_distance_3_object();
  typename K::Construct_supporting_plane_3 supporting_plane =
    k.construct_supporting_plane_3_object();
  typename K::Construct_projected_point_3 projection =
    k.construct_projected_point_3_object();

  // Distance from origin to bound
  const FT bound_sq_dist = sq_distance(origin, bound);

  // Project origin on triangle supporting plane
  const Point_3 proj = projection(supporting_plane(triangle), origin);

  // If point is projected outside, return bound
  if ( compare_sq_distance(origin, proj, bound_sq_dist) == CGAL::LARGER )
  {
    return bound;
  }

  Point_3 moved_point;
  bool inside = is_inside_triangle_3(proj,triangle,moved_point,k);

  // If proj is inside triangle, return it
  if ( inside )
  {
    return proj;
  }

  // Else return the constructed point (nearest point of triangle from proj)
  // if it is closest to origin than bound
  if ( compare_sq_distance(origin, moved_point, bound_sq_dist)
                                                        == CGAL::LARGER )
  {
    return bound;
  }

  return moved_point;
}

}  // end namespace internal


template <class K>
inline
Point_3<K>
nearest_point_3(const Point_3<K>& origin,
                const Triangle_3<K>& triangle,
                const Point_3<K>& bound)
{
  return internal::nearest_point_3(origin, triangle, bound, K());
}

}  // end namespace CGAL


#endif // NEAREST_POINT_TRIANGLE_3_H_
