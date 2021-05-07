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

#ifndef CGAL_DISTANCE_2_POINT_2_LINE_2_H
#define CGAL_DISTANCE_2_POINT_2_LINE_2_H

#include <CGAL/Rational_traits.h>
#include <CGAL/number_utils.h>
#include <CGAL/tags.h>

#include <CGAL/Line_2.h>
#include <CGAL/Point_2.h>

namespace CGAL {
namespace internal {

template <class K>
typename K::FT
squared_distance(const typename K::Point_2& pt,
                 const typename K::Line_2& line,
                 const K&,
                 const Homogeneous_tag&)
{
  typedef typename K::RT RT;
  typedef typename K::FT FT;

  const RT& a = line.a();
  const RT& b = line.b();
  const RT& w = pt.hw();
  RT n = a*pt.hx() + b*pt.hy() + w*line.c();
  RT d = (CGAL_NTS square(a) + CGAL_NTS square(b)) * CGAL_NTS square(w);

  return Rational_traits<FT>().make_rational(CGAL_NTS square(n), d);
}

template <class K>
typename K::FT
squared_distance(const typename K::Point_2& pt,
                 const typename K::Line_2& line,
                 const K&,
                 const Cartesian_tag&)
{
  typedef typename K::FT FT;

  const FT& a = line.a();
  const FT& b = line.b();
  FT n = a*pt.x() + b*pt.y() + line.c();
  FT d = CGAL_NTS square(a) + CGAL_NTS square(b);

  return CGAL_NTS square(n)/d;
}

template <class K>
typename K::FT
squared_distance(const typename K::Point_2& pt,
                 const typename K::Line_2& line,
                 const K& k)
{
  typedef typename K::Kernel_tag Tag;
  Tag tag;
  return squared_distance(pt, line, k, tag);
}

template <class K>
inline typename K::FT
squared_distance(const typename K::Line_2& line,
                 const typename K::Point_2& pt,
                 const K& k)
{
  return internal::squared_distance(pt, line, k);
}

} // namespace internal

template <class K>
inline typename K::FT
squared_distance(const Point_2<K>& pt,
                 const Line_2<K>& line)
{
  return K().compute_squared_distance_2_object()(pt, line);
}

template <class K>
inline typename K::FT
squared_distance(const Line_2<K>& line,
                 const Point_2<K>& pt)
{
  return K().compute_squared_distance_2_object()(line, pt);
}

} // namespace CGAL

#endif // CGAL_DISTANCE_2_POINT_2_LINE_2_H
