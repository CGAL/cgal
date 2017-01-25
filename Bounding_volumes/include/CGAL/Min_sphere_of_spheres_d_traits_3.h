// Copyright (c) 1997  ETH Zurich (Switzerland).
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
// 
//
// Author(s)     : Kaspar Fischer

#ifndef CGAL_MIN_SPHERE_OF_SPHERES_D_TRAITS_3_H
#define CGAL_MIN_SPHERE_OF_SPHERES_D_TRAITS_3_H

#include <CGAL/license/Bounding_volumes.h>


#include <CGAL/tags.h>

namespace CGAL {

 struct Farthest_first_heuristic;

  template<typename K_,                      // kernel
    typename FT_,                            // number type
    typename UseSqrt_ = Tag_false,           // whether to use square-roots
    typename Algorithm_ = Farthest_first_heuristic> // algorithm to use
  class Min_sphere_of_spheres_d_traits_3 {
  public: // types:
    typedef FT_ FT;
    typedef FT_ Radius;
    typedef typename K_::Point_3 Point;
    typedef std::pair<Point,Radius> Sphere;
    typedef typename K_::Cartesian_const_iterator_3 Cartesian_const_iterator;
    typedef UseSqrt_ Use_square_roots;
    typedef Algorithm_ Algorithm;

  public: // constants:
    static const int D = 3;                  // dimension

  public: // accessors:
    static inline const FT& radius(const Sphere& s) {
      return s.second;
    }

    static inline Cartesian_const_iterator
      center_cartesian_begin(const Sphere& s) {
      return s.first.cartesian_begin();
    }
  };

} // namespace CGAL

#endif // CGAL_MIN_SPHERE_OF_SPHERES_D_TRAITS_3_H
