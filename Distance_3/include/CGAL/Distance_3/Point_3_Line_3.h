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
// Author(s)     : Geert-Jan Giezeman, Andreas Fabri, Mael Rouxel-Labbé

#ifndef CGAL_DISTANCE_3_POINT_3_LINE_3_H
#define CGAL_DISTANCE_3_POINT_3_LINE_3_H

#include <CGAL/Distance_3/internal/squared_distance_utils_3.h>

#include <CGAL/Line_3.h>
#include <CGAL/Point_3.h>

namespace CGAL {
namespace internal {

template <class K>
void
squared_distance_RT(const typename K::Point_3 &pt,
                    const typename K::Line_3 &line,
                    typename K::RT& num,
                    typename K::RT& den,
                    const K& k)
{
  typedef typename K::Vector_3 Vector_3;

  typename K::Construct_vector_3 vector = k.construct_vector_3_object();

  const Vector_3& dir = line.direction().vector();
  const Vector_3 diff = vector(line.point(), pt);

  return internal::squared_distance_to_line_RT(dir, diff, num, den, k);
}

template <class K>
typename K::FT
squared_distance(const typename K::Point_3& pt,
                 const typename K::Line_3& line,
                 const K& k)
{
  typedef typename K::RT RT;
  typedef typename K::FT FT;
  RT num, den;
  squared_distance_RT(pt, line, num, den, k);
  return Rational_traits<FT>().make_rational(num, den);
}

template <class K>
inline
typename K::FT
squared_distance(const typename K::Line_3& line,
                 const typename K::Point_3& pt,
                 const K& k)
{
  return squared_distance(pt, line, k);
}

} // namespace internal

template <class K>
inline
typename K::FT
squared_distance(const Point_3<K>& pt,
                 const Line_3<K>& line)
{
  return K().compute_squared_distance_3_object()(pt, line);
}

template <class K>
inline
typename K::FT
squared_distance(const Line_3<K>& line,
                 const Point_3<K>& pt)
{
  return K().compute_squared_distance_3_object()(line, pt);
}

} // namespace CGAL

#endif // CGAL_DISTANCE_3_POINT_3_LINE_3_H
