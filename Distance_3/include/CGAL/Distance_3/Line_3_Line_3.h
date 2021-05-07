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

#ifndef CGAL_DISTANCE_3_LINE_3_LINE_3_H
#define CGAL_DISTANCE_3_LINE_3_LINE_3_H

#include <CGAL/Distance_3/internal/squared_distance_utils_3.h>

#include <CGAL/Line_3.h>

namespace CGAL {
namespace internal {

template <class K>
typename K::FT
squared_distance(const typename K::Line_3& line1,
                 const typename K::Line_3& line2,
                 const K& k)
{
  typedef typename K::Vector_3 Vector_3;

  typename K::Construct_vector_3 vector = k.construct_vector_3_object();

  const Vector_3 dir1 = line1.direction().vector();
  const Vector_3 dir2 = line2.direction().vector();
  const Vector_3 normal = wcross(dir1, dir2, k);
  const Vector_3 diff = vector(line1.point(), line2.point());

  if (is_null(normal, k))
    return squared_distance_to_line(dir2, diff, k);

  return squared_distance_to_plane(normal, diff, k);
}

} // namespace internal

template <class K>
inline
typename K::FT
squared_distance(const Line_3<K>& line1,
                 const Line_3<K>& line2)
{
  return K().compute_squared_distance_3_object()(line1, line2);
}

} // namespace CGAL

#endif // CGAL_DISTANCE_3_LINE_3_LINE_3_H
