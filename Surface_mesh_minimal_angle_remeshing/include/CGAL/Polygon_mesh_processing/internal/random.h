// Copyright (c) 2019  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Kaimo Hu

#ifndef _RANDOM_H_
#define _RANDOM_H_

class Random
{
public:
static double random_in(const double min,
                        const double max) {
  //srand(time(NULL));
  double v = (double)rand() / (double)RAND_MAX;
  v *= max - min;
  return v + min;
}

static double random_sign() {
  double r = random_in(0.0,1.0);
  return r >= 0.5 ? 1.0 : -1.0;
}

//template<class Vector, class Point, class Ray>
//static Ray random_ray_from(const Point& p) {
//  Vector vec = random_vec<Vector>(1.0);
//  return Ray(source,vec);
//}

template<class Vector>
static Vector random_vec(const double mag = 1.0) {
  double x = random_in(-mag,mag);
  double y = random_in(-mag,mag);
  double z = random_in(-mag,mag);
  return Vector(x,y,z);
}

template<class Point>
static Point random_point(const double cmin,
                          const double cmax) {
  double x = random_in(cmin,cmax);
  double y = random_in(cmin,cmax);
  double z = random_in(cmin,cmax);
  return Point(x,y,z);
}

template <class FT>
FT random_value(const FT min, const FT max) {
	FT range = max - min;
  //srand(time(NULL));
	return min + (((FT)rand()) / ((FT)RAND_MAX)) * range;
}

};

#endif // _RANDOM_H_


