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

#ifndef CGAL_DISTANCE_2_LINE_2_LINE_2_H
#define CGAL_DISTANCE_2_LINE_2_LINE_2_H

#include <CGAL/Distance_2/internal/squared_distance_utils_2.h>

#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Line_2.h>

namespace CGAL {
namespace internal {

template <class K>
inline typename K::FT
squared_distance(const typename K::Line_2& line1,
                 const typename K::Line_2& line2,
                 const K& k)
{
  typedef typename K::FT FT;

  typename K::Compute_squared_distance_2 sq_dist = k.compute_squared_distance_2_object();

  if(internal::parallel(line1, line2, k))
    return sq_dist(line1.point(), line2);
  else
    return FT(0);
}

} // namespace internal

template <class K>
inline typename K::FT
squared_distance(const Line_2<K>& line1,
                 const Line_2<K>& line2)
{
  return K().compute_squared_distance_2_object()(line1, line2);
}

} // namespace CGAL

#endif // CGAL_DISTANCE_2_LINE_2_LINE_2_H
