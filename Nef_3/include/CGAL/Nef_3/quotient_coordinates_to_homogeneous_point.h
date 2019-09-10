// Copyright (c) 1997-2000  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
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
