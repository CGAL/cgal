// Copyright (c) 1998-2004
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
// Author(s)     : Geert-Jan Giezeman


#ifndef CGAL_SQUARED_DISTANCE_2_2_H
#define CGAL_SQUARED_DISTANCE_2_2_H

#include <CGAL/user_classes.h>

#include <CGAL/kernel_assertions.h>
#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/Ray_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/enum.h>
#include <CGAL/wmult.h>
#include <CGAL/squared_distance_utils.h>
#include <CGAL/squared_distance_2_1.h>

namespace CGAL {

namespace internal {

  template <class K>
  void
  distance_index(int &ind1,
                 int &ind2,
                 const typename K::Point_2 &pt,
                 const typename K::Triangle_2 &triangle,
                 const K& k )
  {
    typename K::Left_turn_2 leftturn = k.left_turn_2_object();
    typedef typename K::Point_2 Point_2;
    const Point_2 &vt0 = triangle.vertex(0);
    const Point_2 &vt1 = triangle.vertex(1);
    const Point_2 &vt2 = triangle.vertex(2);
    if (leftturn(vt0, vt1, vt2)) {
      if (leftturn(pt, vt1, vt0)) {
        if (!is_acute_angle(vt0, vt1, pt, k)) {
          if (leftturn(pt, vt2, vt1)) {
            if (!is_acute_angle(vt1, vt2, pt, k)) {
              ind1 = 2; ind2 = -1;
              return;
            }
            if (!is_acute_angle(vt2, vt1, pt, k)) {
              ind1 = 1; ind2 = -1;
              return;
            }
            ind1 = 1; ind2 = 2;
            return;
          }
          ind1 = 1; ind2 = -1;
          return;
        }
        if (!is_acute_angle(vt1, vt0, pt, k)) {
          if (leftturn(pt, vt0, vt2)) {
            if (!is_acute_angle(vt0, vt2, pt, k)) {
              ind1 = 2; ind2 = -1;
              return;
            }
            if (!is_acute_angle(vt2, vt0, pt, k)) {
              ind1 = 0; ind2 = -1;
              return;
            }
            ind1 = 2; ind2 = 0;
            return;
          }
          ind1 = 0; ind2 = -1;
          return;
        }
        ind1 = 0; ind2 = 1;
        return;
      } else {
        if (leftturn(pt, vt2, vt1)) {
          if (!is_acute_angle(vt1, vt2, pt, k)) {
            if (leftturn(pt, vt0, vt2)) {
              if (!is_acute_angle(vt0, vt2, pt, k)) {
                ind1 = 2; ind2 = -1;
                return;
              }
              if (!is_acute_angle(vt2, vt0, pt, k)) {
                ind1 = 0; ind2 = -1;
                return;
              }
              ind1 = 2; ind2 = 0;
              return;
            }
            ind1 = 0; ind2 = -1;
            return;
          }
          if (!is_acute_angle(vt2, vt1, pt, k)) {
            ind1 = 1; ind2 = -1;
            return;
          }
          ind1 = 1; ind2 = 2;
          return;
        } else {
          if (leftturn(pt, vt0, vt2)) {
            if (!is_acute_angle(vt2, vt0, pt, k)) {
              ind1 = 0; ind2 = -1;
              return;
            }
            if (!is_acute_angle(vt0, vt2, pt, k)) {
              ind1 = 2; ind2 = -1;
              return;
            }
            ind1 = 2; ind2 = 0;
            return;
          } else {
            ind1 = -1; ind2 = -1; // point inside or on boundary.
            return;
          }
        }
      }
    } else {
      if (leftturn(pt, vt2, vt0)) {
        if (!is_acute_angle(vt0, vt2, pt, k)) {
          if (leftturn(pt, vt1, vt2)) {
            if (!is_acute_angle(vt2, vt1, pt, k)) {
              ind1 = 1; ind2 = -1;
              return;
            }
            if (!is_acute_angle(vt1, vt2, pt, k)) {
              ind1 = 2; ind2 = -1;
              return;
            }
            ind1 = 2; ind2 = 1;
            return;
          }
          ind1 = 2; ind2 = -1;
          return;
        }
        if (!is_acute_angle(vt2, vt0, pt, k)) {
          if (leftturn(pt, vt0, vt1)) {
            if (!is_acute_angle(vt0, vt1, pt, k)) {
              ind1 = 1; ind2 = -1;
              return;
            }
            if (!is_acute_angle(vt1, vt0, pt, k)) {
              ind1 = 0; ind2 = -1;
              return;
            }
            ind1 = 1; ind2 = 0;
            return;
          }
          ind1 = 0; ind2 = -1;
          return;
        }
        ind1 = 0; ind2 = 2;
        return;
      } else {
        if (leftturn(pt, vt1, vt2)) {
          if (!is_acute_angle(vt2, vt1, pt, k)) {
            if (leftturn(pt, vt0, vt1)) {
              if (!is_acute_angle(vt0, vt1, pt, k)) {
                ind1 = 1; ind2 = -1;
                return;
              }
              if (!is_acute_angle(vt1, vt0, pt, k)) {
                ind1 = 0; ind2 = -1;
                return;
              }
              ind1 = 1; ind2 = 0;
              return;
            }
            ind1 = 0; ind2 = -1;
            return;
          }
          if (!is_acute_angle(vt1, vt2, pt, k)) {
            ind1 = 2; ind2 = -1;
            return;
          }
          ind1 = 2; ind2 = 1;
          return;
        } else {
          if (leftturn(pt, vt0, vt1)) {
            if (!is_acute_angle(vt1, vt0, pt, k)) {
              ind1 = 0; ind2 = -1;
              return;
            }
            if (!is_acute_angle(vt0, vt1, pt, k)) {
              ind1 = 1; ind2 = -1;
              return;
            }
            ind1 = 1; ind2 = 0;
            return;
          } else {
            ind1 = -1; ind2 = -1; // point inside or on boundary.
            return;
          }
        }
      }
    }
  }




  template <class K>
  typename K::FT
  squared_distance_indexed(const typename K::Point_2 &pt,
                           const typename K::Triangle_2 &triangle,
                           int ind1, int ind2,
                           const K& k)
  {
    typedef typename K::FT      FT;
    typedef typename K::Line_2  Line_2;
    if (ind1 == -1)
      return FT(0);
    if (ind2 == -1)
      return internal::squared_distance(pt, triangle.vertex(ind1), k);
    return internal::squared_distance(pt,
                 Line_2(triangle.vertex(ind1), triangle.vertex(ind2)),
                                   k);
  }



  template <class K>
  typename K::FT
  squared_distance(const typename K::Point_2 &pt,
                   const typename K::Triangle_2 &triangle,
                   const K& k)
  {
    int ind1,ind2;
    distance_index<K>(ind1, ind2, pt, triangle, k);
    return squared_distance_indexed(pt, triangle, ind1, ind2, k);
  }


  template <class K>
  inline typename K::FT
  squared_distance(const typename K::Triangle_2 & triangle,
                   const typename K::Point_2 & pt,
                   const K& k)
  {
    return internal::squared_distance(pt, triangle, k);
  }


  template <class K>
  typename K::FT
  squared_distance(const typename K::Line_2 &line,
                   const typename K::Triangle_2 &triangle,
                   const K& k)
  {
    typedef typename K::FT FT;
    Oriented_side side0;
    side0 = line.oriented_side(triangle.vertex(0));
    if (line.oriented_side(triangle.vertex(1)) != side0)
      return FT(0);
    if (line.oriented_side(triangle.vertex(2)) != side0)
      return FT(0);
    FT mindist, dist;
    int i;
    mindist = internal::squared_distance(triangle.vertex(0),line,k);
    for (i=1; i<3; i++) {
      dist = internal::squared_distance(triangle.vertex(i),line,k);
      if (dist < mindist)
        mindist = dist;
    }
    return mindist;
  }


  template <class K>
  inline typename K::FT
  squared_distance(const typename K::Triangle_2 & triangle,
                   const typename K::Line_2 & line,
                   const K& k)
  {
    return internal::squared_distance(line, triangle, k);
  }


  template <class K>
  typename K::FT
  squared_distance(const typename K::Ray_2 &ray,
                   const typename K::Triangle_2 &triangle,
                   const K& k)
  {
    typedef typename K::FT       FT;
    typedef typename K::Point_2  Point_2;
    typedef typename K::Line_2   Line_2;
    int i, ind_tr1, ind_tr2, ind_ray = 0, ind1;
    FT mindist, dist;
    distance_index<K>(ind_tr1, ind_tr2, ray.source(), triangle, k);
    mindist =
      squared_distance_indexed(ray.source(), triangle, ind_tr1, ind_tr2, k);
    for (i=0; i<3; i++) {
      const Point_2& pt = triangle.vertex(i);
      distance_index<K>(ind1, pt, ray, k);
      dist = squared_distance_indexed(pt, ray, ind1, k);
      if (dist < mindist) {
        ind_ray = ind1;
        ind_tr1 = i; ind_tr2 = -1;
        mindist = dist;
      }
    }
    // now check if all vertices are on the right side of the separating line.
    // In case of vertex-vertex smallest distance this is the case.
    if (ind_tr2 == -1 && ind_ray != -1)
      return mindist;
    if (ind_tr2 != -1) {
      // Check if all the segment vertices lie at the same side of
      // the triangle segment.
      const Point_2 &vt1 = triangle.vertex(ind_tr1);
      const Point_2 &vt2 = triangle.vertex(ind_tr2);
      if (clockwise(ray.direction().vector(), vt2-vt1, k)) {
        mindist = FT(0);
      }
    } else {
      // Check if all the triangle vertices lie
      // at the same side of the segment.
      const Line_2 &sl = ray.supporting_line();
      Oriented_side or_s = sl.oriented_side(triangle.vertex(0));
      for (i=1; i<3; i++) {
        if (sl.oriented_side(triangle.vertex(i)) != or_s) {
          mindist = FT(0);
          break;
        }
      }
    }
    return mindist;
  }


  template <class K>
  inline typename K::FT
  squared_distance(const typename K::Triangle_2 & triangle,
                   const typename K::Ray_2 & ray,
                   const K& k)
  {
    return internal::squared_distance(ray, triangle, k);
  }


  template <class K>
  typename K::FT
  squared_distance(const typename K::Segment_2 &seg,
                   const typename K::Triangle_2 &triangle,
                   const K& k)
  {
    typedef typename K::FT       FT;
    typedef typename K::Point_2  Point_2;
    typename K::Orientation_2 orientation;
    int i, ind_tr1 = 0, ind_tr2 = -1, ind_seg = 0, ind1, ind2;
    FT mindist, dist;
    mindist = internal::squared_distance(seg.source(), triangle.vertex(0), k);
    for (i=0; i<2; i++) {
      const Point_2 &pt = seg.vertex(i);
      distance_index<K>(ind1, ind2, pt, triangle, k);
      dist = internal::squared_distance_indexed(pt, triangle, ind1, ind2, k);
      if (dist < mindist) {
        ind_seg = i;
        ind_tr1 = ind1; ind_tr2 = ind2;
        mindist = dist;
      }
    }
    for (i=0; i<3; i++) {
      const Point_2& pt = triangle.vertex(i);
      distance_index<K>(ind1, pt, seg, k);
      dist = internal::squared_distance_indexed(pt, seg, ind1, k);
      if (dist < mindist) {
        ind_seg = ind1;
        ind_tr1 = i; ind_tr2 = -1;
        mindist = dist;
      }
    }
    // now check if all vertices are on the right side of the separating line.
    // In case of vertex-vertex smallest distance this is the case.
    if (ind_tr2 == -1 && ind_seg != -1)
      return mindist;

    if (ind_tr2 != -1) {
      // Check if all the segment vertices lie at the same side of
      // the triangle segment.
      const Point_2 &vt1 = triangle.vertex(ind_tr1);
      const Point_2 &vt2 = triangle.vertex(ind_tr2);
      Orientation or_s = orientation(vt1, vt2, seg.source());
      if (orientation(vt1, vt2, seg.target()) != or_s) {
        mindist = FT(0);
      }
    } else {
      // Check if all the triangle vertices lie
      // at the same side of the segment.
      const Point_2 &vt1 = seg.source();
      const Point_2 &vt2 = seg.target();
      Orientation or_s = orientation(vt1, vt2, triangle.vertex(0));
      for (i=1; i<3; i++) {
        if (orientation(vt1, vt2, triangle.vertex(i)) != or_s) {
          mindist = FT(0);
          break;
        }
      }
    }
    return mindist;
  }


  template <class K>
  inline typename K::FT
  squared_distance(const typename K::Triangle_2 & triangle,
                   const typename K::Segment_2 & seg,
                   const K& k)
  {
    return internal::squared_distance(seg, triangle, k);
  }



  template <class K>
  typename K::FT
  squared_distance(const typename K::Triangle_2 &triangle1,
                   const typename K::Triangle_2 &triangle2,
                   const K& k)
  {
    typedef typename K::FT       FT;
    typedef typename K::Point_2  Point_2;
    typename K::Orientation_2 orientation;
    int i, ind1_1 = 0,ind1_2 = -1, ind2_1 = 0, ind2_2 = -1, ind1, ind2;
    FT mindist, dist;
    mindist =
      internal::squared_distance(triangle1.vertex(0), triangle2.vertex(0), k);
    for (i=0; i<3; i++) {
      const Point_2& pt = triangle1.vertex(i);
      distance_index<K>(ind1, ind2, pt, triangle2, k);
      dist = squared_distance_indexed(pt, triangle2, ind1, ind2, k);
      if (dist < mindist) {
        ind1_1 = i; ind1_2 = -1;
        ind2_1 = ind1; ind2_2 = ind2;
        mindist = dist;
      }
    }
    for (i=0; i<3; i++) {
      const Point_2& pt = triangle2.vertex(i);
      distance_index<K>(ind1, ind2, pt, triangle1, k);
      dist = squared_distance_indexed(pt, triangle1, ind1, ind2, k);
      if (dist < mindist) {
        ind1_1 = ind1; ind1_2 = ind2;
        ind2_1 = i; ind2_2 = -1;
        mindist = dist;
      }
    }
    // now check if all vertices are on the right side of the
    // separating line.
    if (ind1_2 == -1 && ind2_2 == -1)
      return mindist;
    // In case of point-segment closest distance, there is still the
    // possibility of overlapping triangles.  Check if all the
    // vertices lie at the same side of the segment.
    if (ind1_2 != -1) {
      const Point_2 &vt1 = triangle1.vertex(ind1_1);
      const Point_2 &vt2 = triangle1.vertex(ind1_2);
      Orientation or_s = orientation(vt1, vt2, triangle2.vertex(0));
      for (i=1; i<3; i++) {
        if (orientation(vt1, vt2, triangle2.vertex(i)) != or_s) {
          mindist = FT(0);
          break;
        }
      }
    } else {
      const Point_2 &vt1 = triangle2.vertex(ind2_1);
      const Point_2 &vt2 = triangle2.vertex(ind2_2);
      Orientation or_s = orientation(vt1, vt2, triangle1.vertex(0));
      for (i=1; i<3; i++) {
        if (orientation(vt1, vt2, triangle1.vertex(i)) != or_s) {
          mindist = FT(0);
          break;
        }
      }
    }
    return mindist;
  }

} // namespace internal

template <class K>
inline typename K::FT
squared_distance(const Point_2<K> &pt,
                 const Triangle_2<K> &triangle)
{
  return internal::squared_distance(pt, triangle, K());
}

template <class K>
inline typename K::FT
squared_distance(const Triangle_2<K> &triangle,
                 const Point_2<K> &pt)
{
  return internal::squared_distance(pt, triangle, K());
}

template <class K>
inline typename K::FT
squared_distance(const Line_2<K> &line,
                 const Triangle_2<K> &triangle)
{
  return internal::squared_distance(line, triangle, K());
}

template <class K>
inline typename K::FT
squared_distance(const Triangle_2<K> &triangle,
                 const Line_2<K> &line)
{
  return internal::squared_distance(line, triangle, K());
}

template <class K>
inline typename K::FT
squared_distance(const Ray_2<K> &ray,
                 const Triangle_2<K> &triangle)
{
  return internal::squared_distance(ray, triangle, K());
}

template <class K>
inline typename K::FT
squared_distance(const Triangle_2<K> &triangle,
                 const Ray_2<K> &ray)
{
  return internal::squared_distance(ray, triangle, K());
}

template <class K>
inline typename K::FT
squared_distance(const Segment_2<K> &seg,
                 const Triangle_2<K> &triangle)
{
  return internal::squared_distance(seg, triangle, K());
}

template <class K>
inline typename K::FT
squared_distance(const Triangle_2<K> &triangle,
                 const Segment_2<K> &seg)
{
  return internal::squared_distance(seg, triangle, K());
}

template <class K>
inline typename K::FT
squared_distance(const Triangle_2<K> &triangle1,
                 const Triangle_2<K> &triangle2)
{
  return internal::squared_distance(triangle1, triangle2, K());
}

} //namespace CGAL

#endif
