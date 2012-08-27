// Copyright (c) 1997-2001  ETH Zurich (Switzerland).
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
// Author(s)     : Kaspar Fischer <fischerk@inf.ethz.ch>

#ifndef CGAL_APPROX_MIN_ELL_D_TRAITS_D_H
#define CGAL_APPROX_MIN_ELL_D_TRAITS_D_H

namespace CGAL {

  template<typename K_, typename ET_>    // kernel and exact number-type
  struct Approximate_min_ellipsoid_d_traits_d {
    typedef double               FT;     // number type (must be double)
    typedef ET_                  ET;     // number type used for exact
                                         // computations (like i.e. computing
                                         // the exact approximation ratio
                                         // in achieved_epsilon())
    typedef typename K_::Point_d Point;  // point type
    typedef typename K_::Cartesian_const_iterator_d Cartesian_const_iterator;
                                         // iterator over point coordinates
    
    static int dimension(const Point& p)
    {
      return p.dimension();
    }
    
    static Cartesian_const_iterator cartesian_begin(const Point& p)
    {
      return p.cartesian_begin();
    }
  };

} // namespace CGAL

#endif // CGAL_APPROX_MIN_ELL_D_TRAITS_D_H
