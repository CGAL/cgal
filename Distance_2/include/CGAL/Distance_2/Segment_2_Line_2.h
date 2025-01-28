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

#ifndef CGAL_DISTANCE_2_SEGMENT_2_LINE_2_H
#define CGAL_DISTANCE_2_SEGMENT_2_LINE_2_H

#include <CGAL/Distance_2/internal/squared_distance_utils_2.h>
#include <CGAL/Distance_2/Point_2_Line_2.h>
#include <CGAL/Distance_2/Segment_2_Ray_2.h>

#include <CGAL/number_utils.h>
#include <CGAL/Rational_traits.h>

#include <CGAL/Line_2.h>
#include <CGAL/Segment_2.h>

namespace CGAL {
namespace internal {

template <class K>
typename K::FT
_sqd_to_line(const typename K::Vector_2& diff,
             const typename K::RT& wcross,
             const typename K::Vector_2& dir)
{
  typedef typename K::RT RT;
  typedef typename K::FT FT;

  RT numerator = CGAL_NTS square(wcross);
  RT denominator = wmult((K*)0, RT(wdot(dir,dir, K())), diff.hw(), diff.hw());
  return Rational_traits<FT>().make_rational(numerator, denominator);
}

template <class K>
typename K::FT
squared_distance(const typename K::Segment_2& seg,
                 const typename K::Line_2& line,
                 const K& k)
{
  typedef typename K::RT RT;
  typedef typename K::FT FT;
  typedef typename K::Vector_2 Vector_2;
  typedef typename K::Point_2  Point_2;

  typename K::Construct_vector_2 vector = k.construct_vector_2_object();
  typename K::Compute_squared_distance_2 sq_dist = k.compute_squared_distance_2_object();

  const Vector_2 linedir = line.direction().vector();
  const Point_2& linepoint = line.point();
  const Vector_2 startvec = vector(linepoint, seg.source());
  const Vector_2 endvec = vector(linepoint, seg.target());

  if(seg.source() == seg.target())
    return sq_dist(seg.source(), line);

  bool crossing1;
  const RT c1s = wcross(linedir, startvec, k);
  const RT c1e = wcross(linedir, endvec, k);
  if(c1s < RT(0))
  {
    crossing1 = (c1e >= RT(0));
  }
  else
  {
    if(c1e <= RT(0))
      crossing1 = true;
     else
      crossing1 = (c1s == RT(0));
  }

  if(crossing1)
  {
    return FT(0);
  }
  else
  {
    const RT dm = _distance_measure_sub<K>(c1s, c1e, startvec, endvec);
    if(dm <= RT(0))
      return _sqd_to_line<K>(startvec, c1s, linedir);
     else
      return _sqd_to_line<K>(endvec, c1e, linedir);
  }
}

template <class K>
inline typename K::FT
squared_distance(const typename K::Line_2& line,
                 const typename K::Segment_2& seg,
                 const K& k)
{
  return internal::squared_distance(seg, line, k);
}

} // namespace internal

template <class K>
inline typename K::FT
squared_distance(const Segment_2<K>& seg,
                 const Line_2<K>& line)
{
  return K().compute_squared_distance_2_object()(seg, line);
}

template <class K>
inline typename K::FT
squared_distance(const Line_2<K>& line,
                 const Segment_2<K>& seg)
{
  return K().compute_squared_distance_2_object()(line, seg);
}

} // namespace CGAL

#endif // CGAL_DISTANCE_2_SEGMENT_2_LINE_2_H
