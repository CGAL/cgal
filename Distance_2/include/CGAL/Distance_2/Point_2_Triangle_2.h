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
//                 Michel Hoffmann <hoffmann@inf.ethz.ch>
//                 Andreas Fabri <Andreas.Fabri@geometryfactory.com>

#ifndef CGAL_DISTANCE_2_POINT_2_TRIANGLE_2_H
#define CGAL_DISTANCE_2_POINT_2_TRIANGLE_2_H

#include <CGAL/Distance_2/internal/squared_distance_utils_2.h>

#include <CGAL/Point_2.h>
#include <CGAL/Triangle_2.h>

namespace CGAL {
namespace internal {

template <class K>
void
distance_index(int &ind1,
               int &ind2,
               const typename K::Point_2 &pt,
               const typename K::Triangle_2 &triangle,
               const K& k)
{
  typedef typename K::Point_2 Point_2;

  typename K::Left_turn_2 leftturn = k.left_turn_2_object();

  const Point_2& vt0 = triangle.vertex(0);
  const Point_2& vt1 = triangle.vertex(1);
  const Point_2& vt2 = triangle.vertex(2);

  if(leftturn(vt0, vt1, vt2)) {
    if(leftturn(pt, vt1, vt0)) {
      if(!is_acute_angle(vt0, vt1, pt, k)) {
        if(leftturn(pt, vt2, vt1)) {
          if(!is_acute_angle(vt1, vt2, pt, k)) {
            ind1 = 2; ind2 = -1;
            return;
          }
          if(!is_acute_angle(vt2, vt1, pt, k)) {
            ind1 = 1; ind2 = -1;
            return;
          }
          ind1 = 1; ind2 = 2;
          return;
        }
        ind1 = 1; ind2 = -1;
        return;
      }
      if(!is_acute_angle(vt1, vt0, pt, k)) {
        if(leftturn(pt, vt0, vt2)) {
          if(!is_acute_angle(vt0, vt2, pt, k)) {
            ind1 = 2; ind2 = -1;
            return;
          }
          if(!is_acute_angle(vt2, vt0, pt, k)) {
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
      if(leftturn(pt, vt2, vt1)) {
        if(!is_acute_angle(vt1, vt2, pt, k)) {
          if(leftturn(pt, vt0, vt2)) {
            if(!is_acute_angle(vt0, vt2, pt, k)) {
              ind1 = 2; ind2 = -1;
              return;
            }
            if(!is_acute_angle(vt2, vt0, pt, k)) {
              ind1 = 0; ind2 = -1;
              return;
            }
            ind1 = 2; ind2 = 0;
            return;
          }
          ind1 = 0; ind2 = -1;
          return;
        }
        if(!is_acute_angle(vt2, vt1, pt, k)) {
          ind1 = 1; ind2 = -1;
          return;
        }
        ind1 = 1; ind2 = 2;
        return;
      } else {
        if(leftturn(pt, vt0, vt2)) {
          if(!is_acute_angle(vt2, vt0, pt, k)) {
            ind1 = 0; ind2 = -1;
            return;
          }
          if(!is_acute_angle(vt0, vt2, pt, k)) {
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
    if(leftturn(pt, vt2, vt0)) {
      if(!is_acute_angle(vt0, vt2, pt, k)) {
        if(leftturn(pt, vt1, vt2)) {
          if(!is_acute_angle(vt2, vt1, pt, k)) {
            ind1 = 1; ind2 = -1;
            return;
          }
          if(!is_acute_angle(vt1, vt2, pt, k)) {
            ind1 = 2; ind2 = -1;
            return;
          }
          ind1 = 2; ind2 = 1;
          return;
        }
        ind1 = 2; ind2 = -1;
        return;
      }
      if(!is_acute_angle(vt2, vt0, pt, k)) {
        if(leftturn(pt, vt0, vt1)) {
          if(!is_acute_angle(vt0, vt1, pt, k)) {
            ind1 = 1; ind2 = -1;
            return;
          }
          if(!is_acute_angle(vt1, vt0, pt, k)) {
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
      if(leftturn(pt, vt1, vt2)) {
        if(!is_acute_angle(vt2, vt1, pt, k)) {
          if(leftturn(pt, vt0, vt1)) {
            if(!is_acute_angle(vt0, vt1, pt, k)) {
              ind1 = 1; ind2 = -1;
              return;
            }
            if(!is_acute_angle(vt1, vt0, pt, k)) {
              ind1 = 0; ind2 = -1;
              return;
            }
            ind1 = 1; ind2 = 0;
            return;
          }
          ind1 = 0; ind2 = -1;
          return;
        }
        if(!is_acute_angle(vt1, vt2, pt, k)) {
          ind1 = 2; ind2 = -1;
          return;
        }
        ind1 = 2; ind2 = 1;
        return;
      } else {
        if(leftturn(pt, vt0, vt1)) {
          if(!is_acute_angle(vt1, vt0, pt, k)) {
            ind1 = 0; ind2 = -1;
            return;
          }
          if(!is_acute_angle(vt0, vt1, pt, k)) {
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
squared_distance_indexed(const typename K::Point_2& pt,
                         const typename K::Triangle_2& triangle,
                         int ind1, int ind2,
                         const K& k)
{
  typedef typename K::FT      FT;
  typedef typename K::Line_2  Line_2;

  typename K::Compute_squared_distance_2 sq_dist = k.compute_squared_distance_2_object();

  if(ind1 == -1)
    return FT(0);

  if(ind2 == -1)
    return sq_dist(pt, triangle.vertex(ind1));

  return sq_dist(pt, Line_2{triangle.vertex(ind1), triangle.vertex(ind2)});
}

template <class K>
typename K::FT
squared_distance(const typename K::Point_2& pt,
                 const typename K::Triangle_2& triangle,
                 const K& k)
{
  int ind1, ind2;
  distance_index<K>(ind1, ind2, pt, triangle, k);
  return squared_distance_indexed(pt, triangle, ind1, ind2, k);
}

template <class K>
inline typename K::FT
squared_distance(const typename K::Triangle_2& triangle,
                 const typename K::Point_2& pt,
                 const K& k)
{
  return internal::squared_distance(pt, triangle, k);
}

} // namespace internal

template <class K>
inline typename K::FT
squared_distance(const Point_2<K>& pt,
                 const Triangle_2<K>& triangle)
{
  return K().compute_squared_distance_2_object()(pt, triangle);
}

template <class K>
inline typename K::FT
squared_distance(const Triangle_2<K>& triangle,
                 const Point_2<K>& pt)
{
  return K().compute_squared_distance_2_object()(triangle, pt);
}

} // namespace CGAL

#endif // CGAL_DISTANCE_2_POINT_2_TRIANGLE_2_H
