// Copyright (c) 1997-2000  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Miguel Granados <granados@mpi-sb.mpg.de>

#ifndef CGAL_NEF_QUOTIENT_COORDINATES_TO_HOMOGENEOUS_POINT_H
#define CGAL_NEF_QUOTIENT_COORDINATES_TO_HOMOGENEOUS_POINT_H

#include <CGAL/license/Nef_3.h>


namespace CGAL {

template <typename Homogeneous>
typename Homogeneous::Point_3
quotient_coordinates_to_homogeneous_point(
                                  typename Homogeneous::FT x,
                                  typename Homogeneous::FT y,
                                  typename Homogeneous::FT z) {
  typedef typename Homogeneous::Point_3 Point_3;
  if( (x.denominator() == y.denominator()) &&
      (x.denominator() == z.denominator())) {
    Point_3 p( x.numerator(),
               y.numerator(),
               z.numerator(),
               x.denominator());
    return normalized(p);
  }
  else {
    Point_3 p( x.numerator()   * y.denominator() * z.denominator(),
               x.denominator() * y.numerator()   * z.denominator(),
               x.denominator() * y.denominator() * z.numerator(),
               x.denominator() * y.denominator() * z.denominator());
    return normalized(p);
  }
}

} //namespace CGAL

#endif // CGAL_NEF_QUOTIENT_COORDINATES_TO_HOMOGENEOUS_POINT_H
