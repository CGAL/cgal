// Copyright (c) 1997  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Kaspar Fischer

#ifndef CGAL_MIN_SPHERE_OF_SPHERES_D_TRAITS_2_H
#define CGAL_MIN_SPHERE_OF_SPHERES_D_TRAITS_2_H

#include <CGAL/Weighted_point.h>

namespace CGAL {

  template<typename K_,                      // kernel
    typename FT_,                            // number type
    typename UseSqrt_ = Tag_false,           // whether to use square-roots
    typename Algorithm_ = Default_algorithm> // algorithm to use
  class Min_sphere_of_spheres_d_traits_2 {
  public: // types:
    typedef FT_ FT;
    typedef FT_ Weight;
    typedef typename K_::Point_2 Point;
    typedef CGAL::Weighted_point<Point,Weight> Sphere;
    typedef typename Point::Cartesian_const_iterator Cartesian_const_iterator;
    typedef UseSqrt_ Use_square_roots;
    typedef Algorithm_ Algorithm;

  public: // constants:
    static const int D = 2;                  // dimension

  public: // accessors:
    static inline const FT radius(const Sphere& s) {
      return s.weight();
    }

    static inline Cartesian_const_iterator
      center_cartesian_begin(const Sphere& s) {
      return s.cartesian_begin();
    }
  };

} // namespace CGAL

#endif // CGAL_MIN_SPHERE_OF_SPHERES_D_TRAITS_2_H
