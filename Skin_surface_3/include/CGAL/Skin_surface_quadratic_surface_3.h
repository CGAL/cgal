// Copyright (c) 2005 Rijksuniversiteit Groningen (Netherlands)
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
// Author(s)     : Nico Kruithof <Nico@cs.rug.nl>

#ifndef CGAL_SKIN_SURFACE_QUADRATIC_SURFACE_H
#define CGAL_SKIN_SURFACE_QUADRATIC_SURFACE_H

#include <CGAL/license/Skin_surface_3.h>

#include <CGAL/Skin_surface_traits_3.h>

namespace CGAL {

template < class SkinSurfaceQuadraticSurfaceTraits_3 >
class Skin_surface_quadratic_surface_3
{
public:
  typedef SkinSurfaceQuadraticSurfaceTraits_3     K;
  typedef typename K::FT                          FT;
  typedef typename K::Point_3                     Point;
  typedef typename K::Vector_3                    Vector;
  typedef typename K::Segment_3                   Segment;
  typedef typename K::Weighted_point_3            Weighted_point;

  Skin_surface_quadratic_surface_3()
    : dim(-1), p(0,0,0), c(0)
  {
    for (int i=0; i<6; i++) Q[i] = 0;
  }

  Skin_surface_quadratic_surface_3(FT Qinput[], Point p, FT c, int d)
    : dim(10+d), p(p), c(c)
  {
    for (int i=0; i<6; i++) Q[i] = Qinput[i];
  }

  Skin_surface_quadratic_surface_3(Weighted_point wp0, FT s)
    : dim(0), p(wp0.point()), c(-s*(1-s)*wp0.weight())
  {
    CGAL_PROFILER(std::string("Constructor : ") +
                  std::string(CGAL_PRETTY_FUNCTION));
    Q[1] = Q[3] = Q[4] = 0;
    Q[0] = Q[2] = Q[5] = (1-s);
  }

  Skin_surface_quadratic_surface_3(Weighted_point wp0,
                                   Weighted_point wp1,
                                   FT s) : dim(1)
  {
    CGAL_PROFILER(std::string("Constructor : ") +
                  std::string(CGAL_PRETTY_FUNCTION));

    K k;
    p = k.construct_weighted_circumcenter_3_object()(wp0,wp1);
    c = s*(1-s)*k.compute_squared_radius_smallest_orthogonal_sphere_3_object()(wp0,wp1);
    Vector t = k.construct_vector_3_object()(k.construct_point_3_object()(wp1),
                                             k.construct_point_3_object()(wp0));

    FT den = t*t;
    Q[0] = (-  t.x()*t.x()/den + (1-s));

    Q[1] = (-2*t.y()*t.x()/den);
    Q[2] = (-  t.y()*t.y()/den + (1-s));

    Q[3] = (-2*t.z()*t.x()/den);
    Q[4] = (-2*t.z()*t.y()/den);
    Q[5] = (-  t.z()*t.z()/den + (1-s));
  }

  Skin_surface_quadratic_surface_3(Weighted_point wp0,
                                   Weighted_point wp1,
                                   Weighted_point wp2,
                                   FT s) : dim(2)
  {
    CGAL_PROFILER(std::string("Constructor : ") +
                  std::string(CGAL_PRETTY_FUNCTION));

    K k;
    p = k.construct_weighted_circumcenter_3_object()(wp0,wp1,wp2);
    c = s*(1-s)*k.compute_squared_radius_smallest_orthogonal_sphere_3_object()(wp0,wp1,wp2);

    Vector t = K().construct_orthogonal_vector_3_object()(k.construct_point_3_object()(wp0),
                                                          k.construct_point_3_object()(wp1),
                                                          k.construct_point_3_object()(wp2));

    FT den = t*t;
    Q[0] = -(-  t.x()*t.x()/den + s);

    Q[1] = -(-2*t.y()*t.x()/den);
    Q[2] = -(-  t.y()*t.y()/den + s);

    Q[3] = -(-2*t.z()*t.x()/den);
    Q[4] = -(-2*t.z()*t.y()/den);
    Q[5] = -(-  t.z()*t.z()/den + s);

  }
  Skin_surface_quadratic_surface_3(Weighted_point wp0,
                                   Weighted_point wp1,
                                   Weighted_point wp2,
                                   Weighted_point wp3,
                                   FT s) : dim(3)
  {
    CGAL_PROFILER(std::string("Constructor : ") +
                  std::string(CGAL_PRETTY_FUNCTION));

    K k;
    p = k.construct_weighted_circumcenter_3_object()(wp0,wp1,wp2,wp3);
    c = s*(1-s)*k.compute_squared_radius_smallest_orthogonal_sphere_3_object()(wp0,wp1,wp2,wp3);
    Q[1] = Q[3] = Q[4] = 0;
    Q[0] = Q[2] = Q[5] = -s;
  }

  template <class Input_point>
  FT value(Input_point const &x) const
  {
    typedef Cartesian_converter<typename Input_point::R, K> Converter;

    FT vx = Converter()(x.x()) - p.x();
    FT vy = Converter()(x.y()) - p.y();
    FT vz = Converter()(x.z()) - p.z();

    return vx*(Q[0]*vx) +
           vy*(Q[1]*vx+Q[2]*vy) +
           vz*(Q[3]*vx+Q[4]*vy+Q[5]*vz) +
           c;
  }
  template <class Input_point>
  Sign sign(Input_point const &x) const {
    return CGAL_NTS sign(value(x));
  }

  template <class Input_point>
  typename Input_point::R::Vector_3 gradient(Input_point const &x)
  {
    // NGHK: We use the kernel trick again: DOCUMENT!!!!
    typedef Cartesian_converter<typename Input_point::R, K> Converter;
    typedef Cartesian_converter<K, typename Input_point::R> Inv_converter;
    Converter conv;
    Inv_converter inv_conv;
    Point xp = conv(x);
    return inv_conv(compute_gradient(xp));
  }

  /// Construct the intersection point with the segment (p0,p1)
  template <class Input_point>
  Input_point to_surface(const Input_point &p0, const Input_point &p1)
  {
    // NGHK: We use the kernel trick again: DOCUMENT!!!!
    typedef Cartesian_converter<typename Input_point::R, K> Converter;
    Converter conv;

    Input_point pp0 = p0, pp1 = p1, mid;
    double sq_d = to_double(squared_distance(pp0,pp1));

    if (value(conv(pp1)) < value(conv(pp0))) {
      std::swap(pp0, pp1);
    }

    while (sq_d > 1.e-10) {
      mid = midpoint(pp0,pp1);
      if (value(conv(mid)) > 0) {
        pp1 = mid;
      } else {
        pp0 = mid;
      }
      sq_d /= 4;
    }
    return midpoint(pp0,pp1);
  }

private:
  //template <class Input_point>
  Vector compute_gradient(Point const &x)
  {
    FT vx = x.x() - p.x();
    FT vy = x.y() - p.y();
    FT vz = x.z() - p.z();

    return Vector(2*Q[0]*vx + Q[1]*vy + Q[3]*vz,
                  Q[1]*vx + 2*Q[2]*vy + Q[4]*vz,
                  Q[3]*vx + Q[4]*vy + 2*Q[5]*vz);
  }

  int dim;
  FT Q[6];
  Point p;
  FT c;
};

} //namespace CGAL

#endif // CGAL_SKIN_SURFACE_QUADRATIC_SURFACE_H
