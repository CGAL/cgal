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

  const Vector_2 &dir1 = seg1.direction().vector();
  const Vector_2 &dir2 = seg2.direction().vector();

  if(same_direction(dir1, dir2, k))
  {
    if(!is_acute_angle(seg1.source(), seg1.target(), seg2.source(), k))
      return internal::squared_distance(seg1.target(), seg2.source(), k);

    if(!is_acute_angle(seg1.target(), seg1.source(), seg2.target(), k))
      return internal::squared_distance(seg1.source(), seg2.target(), k);
  }
  else
  {
    if(!is_acute_angle(seg1.source(), seg1.target(), seg2.target(), k))
      return internal::squared_distance(seg1.target(), seg2.target(), k);

    if(!is_acute_angle(seg1.target(), seg1.source(), seg2.source(), k))
      return internal::squared_distance(seg1.source(), seg2.source(), k);
  }

  return internal::squared_distance(seg2.source(), seg1.supporting_line(), k);
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

  bool crossing1, crossing2;
  RT c1s, c1e, c2s, c2e;

  if(seg1.source() == seg1.target())
    return internal::squared_distance(seg1.source(), seg2, k);

  if(seg2.source() == seg2.target())
    return internal::squared_distance(seg2.source(), seg1, k);

  Orientation o1s = orientation(seg2.source(), seg2.target(), seg1.source());
  Orientation o1e = orientation(seg2.source(), seg2.target(), seg1.target());
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

  Orientation o2s = orientation(seg1.source(), seg1.target(), seg2.source());
  Orientation o2e = orientation(seg1.source(), seg1.target(), seg2.target());
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

    c2s = CGAL_NTS abs(wcross(seg1.source(), seg1.target(), seg2.source(), k));
    c2e = CGAL_NTS abs(wcross(seg1.source(), seg1.target(), seg2.target(), k));
    Comparison_result dm = CGAL_NTS compare(c2s,c2e);

    if(dm == SMALLER)
    {
      return internal::squared_distance(seg2.source(), seg1, k);
    }
    else
    {
      if(dm == LARGER)
      {
        return internal::squared_distance(seg2.target(), seg1, k);
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
    c1s = CGAL_NTS abs(wcross(seg2.source(), seg2.target(), seg1.source(), k));
    c1e = CGAL_NTS abs(wcross(seg2.source(), seg2.target(), seg1.target(), k));
    Comparison_result dm = CGAL_NTS compare(c1s,c1e);
    if(crossing2)
    {
      if(dm == SMALLER)
      {
        return internal::squared_distance(seg1.source(), seg2, k);
      }
      else
      {
        if(dm == LARGER)
          return internal::squared_distance(seg1.target(), seg2, k);
        else  // parallel, should not happen (no crossing)
          return internal::squared_distance_parallel(seg1, seg2, k);
      }
    }
    else
    {
      FT min1, min2;

      if(dm == EQUAL)
        return internal::squared_distance_parallel(seg1, seg2, k);

      min1 = (dm == SMALLER) ? internal::squared_distance(seg1.source(), seg2, k):
                               internal::squared_distance(seg1.target(), seg2, k);

      c2s = CGAL_NTS abs(wcross(seg1.source(), seg1.target(), seg2.source(), k));
      c2e = CGAL_NTS abs(wcross(seg1.source(), seg1.target(), seg2.target(), k));
      dm = CGAL_NTS compare(c2s,c2e);

      if(dm == EQUAL)  // should not happen.
        return internal::squared_distance_parallel(seg1, seg2, k);

      min2 = (dm == SMALLER) ? internal::squared_distance(seg2.source(), seg1, k):
                               internal::squared_distance(seg2.target(), seg1, k);

      return (min1 < min2) ? min1 : min2;
    }
  }
}

template <class K>
typename K::FT
squared_distance(const typename K::Segment_2 &seg1,
                 const typename K::Segment_2 &seg2,
                 const K& k,
                 const Homogeneous_tag&)
{
  typedef typename K::RT RT;
  typedef typename K::FT FT;

  bool crossing1, crossing2;
  RT c1s, c1e, c2s, c2e;

  if(seg1.source() == seg1.target())
    return internal::squared_distance(seg1.source(), seg2, k);

  if(seg2.source() == seg2.target())
    return internal::squared_distance(seg2.source(), seg1, k);

  c1s = wcross(seg2.source(), seg2.target(), seg1.source(), k);
  c1e = wcross(seg2.source(), seg2.target(), seg1.target(), k);
  c2s = wcross(seg1.source(), seg1.target(), seg2.source(), k);
  c2e = wcross(seg1.source(), seg1.target(), seg2.target(), k);

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

    RT dm;
    dm = _distance_measure_sub<K>(c2s,c2e, seg2.source(), seg2.target());
    if(dm < RT(0))
    {
      return internal::squared_distance(seg2.source(), seg1, k);
    }
    else
    {
      if(dm > RT(0))
        return internal::squared_distance(seg2.target(), seg1, k);
      else // parallel, should not happen (no crossing)
        return internal::squared_distance_parallel(seg1, seg2, k);
    }
  }
  else
  {
    if(crossing2)
    {
      RT dm = _distance_measure_sub<K>(c1s, c1e,seg1.source(),seg1.target());
      if(dm < RT(0))
      {
        return internal::squared_distance(seg1.source(), seg2, k);
      }
      else
      {
        if(dm > RT(0))
          return internal::squared_distance(seg1.target(), seg2, k);
        else // parallel, should not happen (no crossing)
          return internal::squared_distance_parallel(seg1, seg2, k);
      }
    }
    else
    {
      FT min1, min2;
      RT dm = _distance_measure_sub<K>(c1s, c1e, seg1.source(), seg1.target());
      if(dm == RT(0))
        return internal::squared_distance_parallel(seg1, seg2, k);

      min1 = (dm < RT(0)) ? internal::squared_distance(seg1.source(), seg2, k):
                            internal::squared_distance(seg1.target(), seg2, k);
      dm = _distance_measure_sub<K>(c2s, c2e, seg2.source(), seg2.target());

      if(dm == RT(0))  // should not happen.
        return internal::squared_distance_parallel(seg1, seg2, k);

      min2 = (dm < RT(0)) ? internal::squared_distance(seg2.source(), seg1, k):
                            internal::squared_distance(seg2.target(), seg1, k);

      return (min1 < min2) ? min1 : min2;
    }
  }
}

} // namespace internal

template <class K>
inline typename K::FT
squared_distance(const Segment_2<K>& seg1,
                 const Segment_2<K>& seg2)
{
  typedef typename K::Kernel_tag Tag;
  return internal::squared_distance(seg1, seg2, K(), Tag());
}

} // namespace CGAL

#endif // CGAL_DISTANCE_2_SEGMENT_2_SEGMENT_2_H
