// Copyright (c) 1997-2002  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Seel    <seel@mpi-sb.mpg.de>
//                 Miguel Granados <granados@mpi-sb.mpg.de>
//                 Susan Hert      <hert@mpi-sb.mpg.de>
//                 Lutz Kettner    <kettner@mpi-sb.mpg.de>
#ifndef CGAL_BOUNDED_SIDE_3_H
#define CGAL_BOUNDED_SIDE_3_H

#include <CGAL/license/Nef_3.h>


#include <CGAL/assertions.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Projection_traits_yz_3.h>
#include <CGAL/Projection_traits_xz_3.h>
#include <CGAL/Projection_traits_xy_3.h>

namespace CGAL {

template <class IteratorForward, class R>
Bounded_side bounded_side_3(IteratorForward first,
                            IteratorForward last,
                            const Point_3<R>& point,
                            const Plane_3<R>& plane)
{
  typedef typename CGAL::Projection_traits_yz_3<R> YZ;
  typedef typename CGAL::Projection_traits_xz_3<R> XZ;
  typedef typename CGAL::Projection_traits_xy_3<R> XY;

  CGAL_assertion(!plane.is_degenerate());

  typename R::Non_zero_dimension_3 non_zero_dimension_3;
  int dir = non_zero_dimension_3(plane.orthogonal_vector());

  if(dir == 0) {
    return bounded_side_2(first, last, point, YZ());
  } else if(dir == 1) {
    return bounded_side_2(first, last, point, XZ());
  } else {
    return bounded_side_2(first, last, point, XY());
  }
}

} //namespace CGAL


#endif // CGAL_BOUNDED_SIDE_3_H
