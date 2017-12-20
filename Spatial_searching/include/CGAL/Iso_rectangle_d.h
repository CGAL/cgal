// Copyright (c) 2002  Utrecht University (The Netherlands).
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
// Authors       : Hans Tangelder (<hanst@cs.uu.nl>)

#ifndef CGAL_ISO_RECTANGLE_D_H
#define CGAL_ISO_RECTANGLE_D_H

#include <CGAL/license/Spatial_searching.h>


#include <functional>
#include <algorithm>
#include <new>
#include <CGAL/assertions.h>

namespace CGAL {
  
  template <class R> class Iso_rectangle_d {
  public:
    typedef typename R::FT FT;
    typedef typename R::Point_d Point_d;

  private:

    int dim;
    FT *lower;
    FT *upper;
    
  public:

    Iso_rectangle_d(const Point_d& p, const Point_d& q)
    { CGAL_assertion(p.dimension() == q.dimension());
      dim = p.dimension();
      lower = new FT[dim];
      upper = new FT[dim];
      for (int i = 0; i < dim; ++i) {
	  if (p[i] <= q[i]) {
		lower[i]=p[i]; 
                upper[i]=q[i];
	  }
	  else {
		lower[i]=q[i]; 
                upper[i]=p[i];
	  }
     }	  
    }
  
  // copy constructor
  Iso_rectangle_d(const Iso_rectangle_d& b) : dim(b.dim) {
      lower = new FT[dim];
      upper = new FT[dim];
      for (int i = 0; i < dim; ++i) {
		lower[i]=b.lower[i]; 
                upper[i]=b.upper[i];
      }
  }

  
  bool has_on_bounded_side(const Point_d& p) const
  {
    FT h;
    for (int i = 0; i < dimension(); ++i) {
        h=p[i];
        if ( (h < lower[i]) || (h > upper[i]) ) return 0;
    }
    return 1;
  } 

    inline int dimension() const { return dim;}
    
    inline FT min_coord(int i) const {
      return lower[i];
    }

    inline FT max_coord(int i) const {
      return upper[i];
    }

  }; // end of class

} // namespace CGAL

#endif // CGAL_ISO_RECTANGLE_D_H
