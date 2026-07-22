// Copyright (c) 1998-2021
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_DISTANCE_3_POINT_3_TETRAHEDRON_3_H
#define CGAL_DISTANCE_3_POINT_3_TETRAHEDRON_3_H

#include <CGAL/Distance_3/Point_3_Triangle_3.h>

#include <CGAL/Point_3.h>
#include <CGAL/Tetrahedron_3.h>

namespace CGAL {
namespace internal {

template <class K>
inline
typename K::FT
squared_distance(const typename K::Point_3& pt,
                 const typename K::Tetrahedron_3& tet,
                 const K& k)
{
  typedef typename K::Point_3 Point_3;

  typename K::Construct_vertex_3 vertex = k.construct_vertex_3_object();
  typename K::Orientation_3 orientation = k.orientation_3_object();

  bool on_bounded_side = true;
  const Point_3& t0 = vertex(tet, 0);
  const Point_3& t1 = vertex(tet, 1);
  const Point_3& t2 = vertex(tet, 2);
  const Point_3& t3 = vertex(tet, 3);

  Orientation ori_tet = orientation(t0,t1,t2,t3);

  if(ori_tet==COPLANAR){
    //Degen Tetrahedron
    //Get the minimum of three triangles (no need to test the fourth)
    bool inside = false;
    typename K::FT dmin = squared_distance_to_triangle(pt, t0, t1, t2, k, inside);
    if(inside)
      return dmin;
    typename K::FT d = squared_distance_to_triangle(pt, t0, t1, t3, k, inside);
    if(inside)
      return d;
    dmin=(std::min)(d,dmin);
    d = squared_distance_to_triangle(pt, t0, t2, t3, k, inside);
    return (std::min)(d,dmin);
  }

  bool dmin_initialized = false; //dmin_initialized and !on_bounded_side have always the samed value
  typename K::FT dmin;
  bool inside = false;

  if(orientation(pt, t0,t1,t2) == ori_tet)
  {
    on_bounded_side = false;
    dmin = squared_distance_to_triangle(pt, t0, t1, t2, k, inside);
    dmin_initialized = true;
    if(inside)
      return dmin;
  }

  if(orientation(pt, t0,t3,t1) == ori_tet)
  {
    on_bounded_side = false;
    const typename K::FT d = squared_distance_to_triangle(pt, t0, t3, t1, k, inside);
    if(inside)
      return d;

    if(!dmin_initialized)
    {
      dmin = d;
      dmin_initialized = true;
    }
    else
    {
      dmin = (std::min)(d, dmin);
    }
  }

  if(orientation(pt, t1,t3,t2) == ori_tet)
  {
    on_bounded_side = false;
    const typename K::FT d = squared_distance_to_triangle(pt, t1, t3, t2, k, inside);
    if(inside)
      return d;

    if(!dmin_initialized)
    {
      dmin = d;
      dmin_initialized = true;
    }
    else
    {
      dmin = (std::min)(d, dmin);
    }
  }

  if(orientation(pt, t2,t3,t0) == ori_tet)
  {
    on_bounded_side = false;
    const typename K::FT d = squared_distance_to_triangle(pt, t2, t3, t0, k, inside);
    if(inside)
      return d;

    if(!dmin_initialized)
    {
      dmin = d;
      dmin_initialized = true;
    }
    else
    {
      dmin = (std::min)(d, dmin);
    }
  }

  if(on_bounded_side)
    return typename K::FT(0);

  return dmin;
}

template <class K>
inline
typename K::FT
squared_distance(const typename K::Tetrahedron_3& tet,
                 const typename K::Point_3& pt,
                 const K& k)
{
  return squared_distance(pt, tet, k);
}

template <class K>
typename K::Comparison_result
compare_squared_distance(const typename K::Point_3& pt,
                         const typename K::Tetrahedron_3& tet,
                         const K& k,
                         const typename K::FT& d2)
{
  typedef typename K::Point_3 Point_3;

  typename K::Construct_vertex_3 vertex = k.construct_vertex_3_object();
  typename K::Orientation_3 orientation = k.orientation_3_object();

  /* The content of this function is very similar with the one above, the difference is we can exit earlier if
      we found a triangle closer than d or plane farther than d since we do not need the exact distance.
      (there are also early exits in calling functions)  */

  bool on_bounded_side = true;
  bool inside_or_far_to_the_plane = false;
  const Point_3& t0 = vertex(tet, 0);
  const Point_3& t1 = vertex(tet, 1);
  const Point_3& t2 = vertex(tet, 2);
  const Point_3& t3 = vertex(tet, 3);

  Orientation ori_tet = orientation(t0,t1,t2,t3);

  if(ori_tet==COPLANAR){
    //Degen Tetrahedron
    //Get the minimum of three triangles (no need to test the fourth)
    typename K::Comparison_result res = compare_squared_distance_to_triangle(pt, t0, t1, t2, k, d2, inside_or_far_to_the_plane);
    if(inside_or_far_to_the_plane)
      return res;
    typename K::Comparison_result temp_res = compare_squared_distance_to_triangle(pt, t0, t1, t3, k, d2, inside_or_far_to_the_plane);
    if(inside_or_far_to_the_plane)
      return temp_res;
    res=smaller_of(res,temp_res);
    temp_res = compare_squared_distance_to_triangle(pt, t0, t2, t3, k, d2, inside_or_far_to_the_plane);
    return smaller_of(res,temp_res);
  }

  typename K::Comparison_result res=LARGER;
  if(orientation(pt, t0,t1,t2) == ori_tet)
  {
    on_bounded_side = false;
    res = compare_squared_distance_to_triangle(pt, t0, t1, t2, k, d2, inside_or_far_to_the_plane);
    if(inside_or_far_to_the_plane || res==SMALLER)
      return res;
  }

  if(orientation(pt, t0,t3,t1) == ori_tet)
  {
    on_bounded_side = false;
    const typename K::Comparison_result temp_res = compare_squared_distance_to_triangle(pt, t0, t3, t1, k, d2, inside_or_far_to_the_plane);
    if(inside_or_far_to_the_plane || temp_res==SMALLER)
      return temp_res;
    res = smaller_of(res, temp_res);
  }

  if(orientation(pt, t1,t3,t2) == ori_tet)
  {
    on_bounded_side = false;
    const typename K::Comparison_result temp_res = compare_squared_distance_to_triangle(pt, t1, t3, t2, k, d2, inside_or_far_to_the_plane);
    if(inside_or_far_to_the_plane || temp_res==SMALLER)
      return temp_res;
    res = smaller_of(res, temp_res);
  }

  if(orientation(pt, t2,t3,t0) == ori_tet)
  {
    on_bounded_side = false;
    const typename K::Comparison_result temp_res = compare_squared_distance_to_triangle(pt, t2, t3, t0, k, d2, inside_or_far_to_the_plane);
    if(inside_or_far_to_the_plane || temp_res==SMALLER)
      return temp_res;
    res = smaller_of(res, temp_res);
  }

  if(on_bounded_side)
    return ::CGAL::compare(typename K::FT(0),d2);

  return res;
}

template <class K>
inline
typename K::Comparison_result
compare_squared_distance(const typename K::Tetrahedron_3& tet,
                         const typename K::Point_3& pt,
                         const K& k,
                         const typename K::FT& d2)
{
  return compare_squared_distance(pt, tet, k, d2);
}

} // namespace internal

} // namespace CGAL

#endif // CGAL_DISTANCE_3_POINT_3_TETRAHEDRON_3_H
