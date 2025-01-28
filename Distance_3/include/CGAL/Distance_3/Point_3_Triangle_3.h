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
// Author(s)     : Geert-Jan Giezeman, Andreas Fabri

#ifndef CGAL_DISTANCE_3_POINT_3_TRIANGLE_3_H
#define CGAL_DISTANCE_3_POINT_3_TRIANGLE_3_H

#include <CGAL/Distance_3/internal/squared_distance_utils_3.h>
#include <CGAL/Distance_3/Point_3_Segment_3.h>

#include <CGAL/Point_3.h>
#include <CGAL/Triangle_3.h>

namespace CGAL {
namespace internal {

// returns true iff pt is on the negative side of the plane defined by (ep0, ep1) and normal
template <class K>
inline bool
on_left_of_triangle_edge(const typename K::Point_3& pt,
                         const typename K::Vector_3& normal,
                         const typename K::Point_3& ep0,
                         const typename K::Point_3& ep1,
                         const K& k)
{
  typedef typename K::RT RT;
  typedef typename K::Vector_3 Vector_3;

  typename K::Construct_vector_3 vector = k.construct_vector_3_object();

  const Vector_3 edge = vector(ep0, ep1);
  const Vector_3 diff = vector(ep0, pt);

  return (wdot(wcross(edge, normal, k), diff, k) <= RT(0));
}

template <class K>
inline void
squared_distance_to_triangle_RT(const typename K::Point_3& pt,
                                const typename K::Point_3& t0,
                                const typename K::Point_3& t1,
                                const typename K::Point_3& t2,
                                bool& inside,
                                typename K::RT& num,
                                typename K::RT& den,
                                const K& k)
{
  typedef typename K::Vector_3 Vector_3;

  typename K::Construct_vector_3 vector = k.construct_vector_3_object();
  typename K::Construct_segment_3 segment = k.construct_segment_3_object();

  const Vector_3 e1 = vector(t0, t1);
  const Vector_3 oe3 = vector(t0, t2);
  const Vector_3 normal = wcross(e1, oe3, k);

  if(normal == NULL_VECTOR)
  {
    // The case normal==NULL_VECTOR covers the case when the triangle
    // is colinear, or even more degenerate. In that case, we can
    // simply take also the distance to the three segments.
    squared_distance_RT(pt, segment(t2, t0), num, den, k);

    typename K::RT num2, den2;
    squared_distance_RT(pt, segment(t1, t2), num2, den2, k);
    if(compare_quotients(num2,den2, num,den) == SMALLER)
    {
      num = num2;
      den = den2;
    }

    // Should not be needed since at most 2 edges cover the full triangle in the degenerate case,
    // but leaving it for robustness
    squared_distance_RT(pt, segment(t0, t1), num2, den2, k);
    if(compare_quotients(num2,den2, num,den) == SMALLER)
    {
      num = num2;
      den = den2;
    }

    return;
  }

  if(!on_left_of_triangle_edge(pt, normal, t0, t1, k))
  {
    if(!on_left_of_triangle_edge(pt, normal, t1, t2, k))
    {
      // can't be to the right of all three segments
      squared_distance_RT(pt, segment(t0, t1), num, den, k);

      typename K::RT num2, den2;
      squared_distance_RT(pt, segment(t1, t2), num2, den2, k);
      if(compare_quotients(num2,den2, num,den) == SMALLER)
      {
        num = num2;
        den = den2;
      }
    }
    else // on_left_of_triangle_edge(pt, normal, t1, t2, k)
    {
      if(!on_left_of_triangle_edge(pt, normal, t2, t0, k))
      {
        squared_distance_RT(pt, segment(t0, t1), num, den, k);

        typename K::RT num2, den2;
        squared_distance_RT(pt, segment(t2, t0), num2, den2, k);
        if(compare_quotients(num2,den2, num,den) == SMALLER)
        {
          num = num2;
          den = den2;
        }
      }
      else // on_left_of_triangle_edge(pt, normal, t2, t0, k)
      {
        return squared_distance_RT(pt, segment(t0, t1), num, den, k);
      }
    }
  }
  else // on_left_of_triangle_edge(pt, normal, t0, t1, k)
  {
    if(!on_left_of_triangle_edge(pt, normal, t1, t2, k))
    {
      if(!on_left_of_triangle_edge(pt, normal, t2, t0, k))
      {
        squared_distance_RT(pt, segment(t1, t2), num, den, k);

        typename K::RT num2, den2;
        squared_distance_RT(pt, segment(t2, t0), num2, den2, k);
        if(compare_quotients(num2,den2, num,den) == SMALLER)
        {
          num = num2;
          den = den2;
        }
      }
      else // on_left_of_triangle_edge(pt, normal, t2, t0, k)
      {
        return squared_distance_RT(pt, segment(t1, t2), num, den, k);
      }
    }
    else // on_left_of_triangle_edge(pt, normal, t1, t2, k)
    {
      if(!on_left_of_triangle_edge(pt, normal, t2, t0, k))
      {
        return squared_distance_RT(pt, segment(t2, t0), num, den, k);
      }
      else // on_left_of_triangle_edge(pt, normal, t2, t0, k)
      {
        inside = true;
        return squared_distance_to_plane_RT(normal, vector(t0, pt), num, den, k);
      }
    }
  }
}

template <class K>
void
squared_distance_RT(const typename K::Point_3& pt,
                    const typename K::Triangle_3& t,
                    typename K::RT& num,
                    typename K::RT& den,
                    const K& k)
{
  typename K::Construct_vertex_3 vertex;
  bool inside = false;
  squared_distance_to_triangle_RT(pt,
                                  vertex(t, 0),
                                  vertex(t, 1),
                                  vertex(t, 2),
                                  inside,
                                  num,
                                  den,
                                  k);
}

template <class K>
inline typename K::FT
squared_distance_to_triangle(const typename K::Point_3& pt,
                             const typename K::Point_3& t0,
                             const typename K::Point_3& t1,
                             const typename K::Point_3& t2,
                             const K& k,
                             bool& inside)
{
  typedef typename K::Vector_3 Vector_3;

  typename K::Construct_segment_3 segment = k.construct_segment_3_object();
  typename K::Construct_vector_3 vector = k.construct_vector_3_object();
  typename K::Compute_squared_distance_3 sq_dist = k.compute_squared_distance_3_object();

  const Vector_3 e1 = vector(t0, t1);
  const Vector_3 oe3 = vector(t0, t2);
  const Vector_3 normal = wcross(e1, oe3, k);

  if(normal == NULL_VECTOR)
  {
    // The case normal == NULL_VECTOR covers the case when the triangle
    // is colinear, or even more degenerate. In that case, we can
    // simply take also the distance to the three segments.
    //
    // Note that in the degenerate case, at most 2 edges cover the full triangle,
    // and only two distances could be used, but leaving 3 for the case of
    // inexact constructions as it might improve the accuracy.
    typename K::FT d1 = sq_dist(pt, segment(t2, t0));
    typename K::FT d2 = sq_dist(pt, segment(t1, t2));
    typename K::FT d3 = sq_dist(pt, segment(t0, t1));

    return (std::min)((std::min)(d1, d2), d3);
  }

  if(!on_left_of_triangle_edge(pt, normal, t0, t1, k))
  {
    if(!on_left_of_triangle_edge(pt, normal, t1, t2, k))
    {
      // can't be to the right of all three segments
      return (std::min)(sq_dist(pt, segment(t0, t1)), sq_dist(pt, segment(t1, t2)));
    }
    else // on_left_of_triangle_edge(pt, normal, t1, t2, k)
    {
      if(!on_left_of_triangle_edge(pt, normal, t2, t0, k))
      {
        return (std::min)(sq_dist(pt, segment(t0, t1)), sq_dist(pt, segment(t2, t0)));
      }
      else // on_left_of_triangle_edge(pt, normal, t2, t0, k)
      {
        return sq_dist(pt, segment(t0, t1));
      }
    }
  }
  else // on_left_of_triangle_edge(pt, normal, t0, t1, k)
  {
    if(!on_left_of_triangle_edge(pt, normal, t1, t2, k))
    {
      if(!on_left_of_triangle_edge(pt, normal, t2, t0, k))
      {
        return (std::min)(sq_dist(pt, segment(t1, t2)), sq_dist(pt, segment(t2, t0)));
      }
      else // on_left_of_triangle_edge(pt, normal, t2, t0, k)
      {
        return sq_dist(pt, segment(t1, t2));
      }
    }
    else // on_left_of_triangle_edge(pt, normal, t1, t2, k)
    {
      if(!on_left_of_triangle_edge(pt, normal, t2, t0, k))
      {
        return sq_dist(pt, segment(t2, t0));
      }
      else // on_left_of_triangle_edge(pt, normal, t2, t0, k)
      {
        inside = true;
        return squared_distance_to_plane(normal, vector(t0, pt), k);
      }
    }
  }
}

template <class K>
inline typename K::FT
squared_distance(const typename K::Point_3& pt,
                 const typename K::Triangle_3& t,
                 const K& k)
{
  typename K::Construct_vertex_3 vertex = k.construct_vertex_3_object();

  bool unused_inside = false;
  return squared_distance_to_triangle(pt,
                                      vertex(t, 0),
                                      vertex(t, 1),
                                      vertex(t, 2),
                                      k,
                                      unused_inside);
}

template <class K>
inline typename K::FT
squared_distance(const typename K::Triangle_3& t,
                 const typename K::Point_3& pt,
                 const K& k)
{
  return squared_distance(pt, t, k);
}

} // namespace internal

template <class K>
inline
typename K::FT
squared_distance(const Point_3<K>& pt,
                 const Triangle_3<K>& t)
{
  return K().compute_squared_distance_3_object()(pt, t);
}

template <class K>
inline
typename K::FT
squared_distance(const Triangle_3<K>& t,
                 const Point_3<K>& pt)
{
  return K().compute_squared_distance_3_object()(t, pt);
}

} // namespace CGAL

#endif // CGAL_DISTANCE_3_POINT_3_TRIANGLE_3_H
