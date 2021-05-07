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

#ifndef CGAL_DISTANCE_2_SEGMENT_2_SEGMENT_2_H
#define CGAL_DISTANCE_2_SEGMENT_2_SEGMENT_2_H

#include <CGAL/Distance_2/internal/squared_distance_utils_2.h>
#include <CGAL/Distance_2/Point_2_Line_2.h>
#include <CGAL/Distance_2/Point_2_Point_2.h>
#include <CGAL/Distance_2/Point_2_Segment_2.h>

#include <CGAL/number_utils.h>
#include <CGAL/tags.h>

#include <CGAL/Segment_2.h>

namespace CGAL {
namespace internal {

template <class K>
typename K::FT
squared_distance_parallel(const typename K::Segment_2& seg1,
                          const typename K::Segment_2& seg2,
                          const K& k)
{
  typedef typename K::Vector_2 Vector_2;

  typename K::Compute_squared_distance_2 sq_dist = k.compute_squared_distance_2_object();

  const Vector_2 dir1 = seg1.direction().vector();
  const Vector_2 dir2 = seg2.direction().vector();

  if(same_direction(dir1, dir2, k))
  {
    if(!is_acute_angle(seg1.source(), seg1.target(), seg2.source(), k))
      return sq_dist(seg1.target(), seg2.source());

    if(!is_acute_angle(seg1.target(), seg1.source(), seg2.target(), k))
      return sq_dist(seg1.source(), seg2.target());
  }
  else
  {
    if(!is_acute_angle(seg1.source(), seg1.target(), seg2.target(), k))
      return sq_dist(seg1.target(), seg2.target());

    if(!is_acute_angle(seg1.target(), seg1.source(), seg2.source(), k))
      return sq_dist(seg1.source(), seg2.source());
  }

  return sq_dist(seg2.source(), seg1.supporting_line());
}

template <class K>
inline typename K::RT
_distance_measure_sub(const typename K::RT& startwcross,
                      const typename K::RT& endwcross,
                      const typename K::Point_2& start,
                      const typename K::Point_2& end)
{
  return CGAL_NTS abs(wmult((K*)0, startwcross, end.hw())) -
           CGAL_NTS abs(wmult((K*)0, endwcross, start.hw()));
}

template <class K>
typename K::FT
squared_distance(const typename K::Segment_2& seg1,
                 const typename K::Segment_2& seg2,
                 const K& k,
                 const Cartesian_tag&)
{
  typedef typename K::RT RT;
  typedef typename K::FT FT;

  typename K::Orientation_2 orientation = k.orientation_2_object();
  typename K::Compute_squared_distance_2 sq_dist = k.compute_squared_distance_2_object();

  if(seg1.source() == seg1.target())
    return sq_dist(seg1.source(), seg2);

  if(seg2.source() == seg2.target())
    return sq_dist(seg2.source(), seg1);

  const Orientation o1s = orientation(seg2.source(), seg2.target(), seg1.source());
  const Orientation o1e = orientation(seg2.source(), seg2.target(), seg1.target());

  bool crossing1, crossing2;
  if(o1s == RIGHT_TURN)
  {
    crossing1 = (o1e != RIGHT_TURN);
  }
  else
  {
    if(o1e != LEFT_TURN)
    {
      if(o1s == COLLINEAR && o1e == COLLINEAR)
        return internal::squared_distance_parallel(seg1, seg2, k);

      crossing1 = true;
    }
    else
    {
      crossing1 = (o1s == COLLINEAR);
    }
  }

  const Orientation o2s = orientation(seg1.source(), seg1.target(), seg2.source());
  const Orientation o2e = orientation(seg1.source(), seg1.target(), seg2.target());
  if(o2s == RIGHT_TURN)
  {
    crossing2 = (o2e != RIGHT_TURN);
  }
  else
  {
    if(o2e != LEFT_TURN)
    {
      if(o2s == COLLINEAR && o2e == COLLINEAR)
        return internal::squared_distance_parallel(seg1, seg2, k);

      crossing2 = true;
    }
    else
    {
      crossing2 = (o2s == COLLINEAR);
    }
  }

  if(crossing1)
  {
    if(crossing2)
      return (FT)0;

    const RT c2s = CGAL_NTS abs(wcross(seg1.source(), seg1.target(), seg2.source(), k));
    const RT c2e = CGAL_NTS abs(wcross(seg1.source(), seg1.target(), seg2.target(), k));
    Comparison_result dm = CGAL_NTS compare(c2s, c2e);

    if(dm == SMALLER)
    {
      return sq_dist(seg2.source(), seg1);
    }
    else
    {
      if(dm == LARGER)
      {
        return sq_dist(seg2.target(), seg1);
      }
      else
      {
        // parallel, should not happen (no crossing)
        return internal::squared_distance_parallel(seg1, seg2, k);
      }
    }
  }
  else
  {
    const RT c1s = CGAL_NTS abs(wcross(seg2.source(), seg2.target(), seg1.source(), k));
    const RT c1e = CGAL_NTS abs(wcross(seg2.source(), seg2.target(), seg1.target(), k));
    Comparison_result dm = CGAL_NTS compare(c1s, c1e);
    if(crossing2)
    {
      if(dm == SMALLER)
      {
        return sq_dist(seg1.source(), seg2);
      }
      else
      {
        if(dm == LARGER)
          return sq_dist(seg1.target(), seg2);
        else  // parallel, should not happen (no crossing)
          return internal::squared_distance_parallel(seg1, seg2, k);
      }
    }
    else
    {
      if(dm == EQUAL)
        return internal::squared_distance_parallel(seg1, seg2, k);

      FT min1 = (dm == SMALLER) ? sq_dist(seg1.source(), seg2):
                                  sq_dist(seg1.target(), seg2);

      const RT c2s = CGAL_NTS abs(wcross(seg1.source(), seg1.target(), seg2.source(), k));
      const RT c2e = CGAL_NTS abs(wcross(seg1.source(), seg1.target(), seg2.target(), k));

      dm = CGAL_NTS compare(c2s,c2e);
      if(dm == EQUAL)  // should not happen.
        return internal::squared_distance_parallel(seg1, seg2, k);

      FT min2 = (dm == SMALLER) ? sq_dist(seg2.source(), seg1):
                                  sq_dist(seg2.target(), seg1);

      return (min1 < min2) ? min1 : min2;
    }
  }
}

template <class K>
typename K::FT
squared_distance(const typename K::Segment_2& seg1,
                 const typename K::Segment_2& seg2,
                 const K& k,
                 const Homogeneous_tag&)
{
  typedef typename K::RT RT;
  typedef typename K::FT FT;

  typename K::Compute_squared_distance_2 sq_dist = k.compute_squared_distance_2_object();

  if(seg1.source() == seg1.target())
    return sq_dist(seg1.source(), seg2);

  if(seg2.source() == seg2.target())
    return sq_dist(seg2.source(), seg1);

  const RT c1s = wcross(seg2.source(), seg2.target(), seg1.source(), k);
  const RT c1e = wcross(seg2.source(), seg2.target(), seg1.target(), k);
  const RT c2s = wcross(seg1.source(), seg1.target(), seg2.source(), k);
  const RT c2e = wcross(seg1.source(), seg1.target(), seg2.target(), k);

  bool crossing1, crossing2;

  if(c1s < RT(0))
  {
    crossing1 = (c1e >= RT(0));
  }
  else
  {
    if(c1e <= RT(0))
    {
      if(c1s == RT(0) && c1e == RT(0))
        return internal::squared_distance_parallel(seg1, seg2, k);

      crossing1 = true;
    }
    else
    {
      crossing1 = (c1s == RT(0));
    }
  }

  if(c2s < RT(0))
  {
    crossing2 = (c2e >= RT(0));
  }
  else
  {
    if(c2e <= RT(0))
    {
      if(c2s == RT(0) && c2e == RT(0))
        return internal::squared_distance_parallel(seg1, seg2, k);

      crossing2 = true;
    }
    else
    {
      crossing2 = (c2s == RT(0));
    }
  }

  if(crossing1)
  {
    if(crossing2)
      return (FT)0;

    const RT dm = _distance_measure_sub<K>(c2s,c2e, seg2.source(), seg2.target());
    if(dm < RT(0))
    {
      return sq_dist(seg2.source(), seg1);
    }
    else
    {
      if(dm > RT(0))
        return sq_dist(seg2.target(), seg1);
      else // parallel, should not happen (no crossing)
        return internal::squared_distance_parallel(seg1, seg2, k);
    }
  }
  else
  {
    if(crossing2)
    {
      const RT dm = _distance_measure_sub<K>(c1s, c1e,seg1.source(),seg1.target());
      if(dm < RT(0))
      {
        return sq_dist(seg1.source(), seg2);
      }
      else
      {
        if(dm > RT(0))
          return sq_dist(seg1.target(), seg2);
        else // parallel, should not happen (no crossing)
          return internal::squared_distance_parallel(seg1, seg2, k);
      }
    }
    else
    {
      RT dm = _distance_measure_sub<K>(c1s, c1e, seg1.source(), seg1.target());
      if(dm == RT(0))
        return internal::squared_distance_parallel(seg1, seg2, k);

      FT min1 = (dm < RT(0)) ? sq_dist(seg1.source(), seg2):
                               sq_dist(seg1.target(), seg2);

      dm = _distance_measure_sub<K>(c2s, c2e, seg2.source(), seg2.target());
      if(dm == RT(0))  // should not happen.
        return internal::squared_distance_parallel(seg1, seg2, k);

      FT min2 = (dm < RT(0)) ? sq_dist(seg2.source(), seg1):
                               sq_dist(seg2.target(), seg1);

      return (min1 < min2) ? min1 : min2;
    }
  }
}

template <class K>
inline typename K::FT
squared_distance(const Segment_2<K>& seg1,
                 const Segment_2<K>& seg2,
                 const K& k)
{
  typedef typename K::Kernel_tag Tag;
  return squared_distance(seg1, seg2, k, Tag());
}

} // namespace internal

template <class K>
inline typename K::FT
squared_distance(const Segment_2<K>& seg1,
                 const Segment_2<K>& seg2)
{
  return K().compute_squared_distance_2_object()(seg1, seg2);
}

} // namespace CGAL

#endif // CGAL_DISTANCE_2_SEGMENT_2_SEGMENT_2_H
