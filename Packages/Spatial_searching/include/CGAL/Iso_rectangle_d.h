// ======================================================================
// Copyright (c) 2002 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Iso_rectangle_d.h
// package       : ASPAS
// revision      : 1.4 
// revision_date : 2003/02/01 
// authors       : Hans Tangelder (<hanst@cs.uu.nl>)
// maintainer    : Hans Tangelder (<hanst@cs.uu.nl>)
// coordinator   : Utrecht University
//
// ======================================================================

#ifndef CGAL_ISO_RECTANGLE_D_H
#define CGAL_ISO_RECTANGLE_D_H
#include <functional>
#include <algorithm>
#include <new>

namespace CGAL {
  
  template <class R> class Iso_rectangle_d {
  public:
    typedef typename R::FT NT;
    typedef typename R::Point_d Point_d;

  private:

    int dim;
    NT *lower;
    NT *upper;
    
  public:

    Iso_rectangle_d(const Point_d& p, const Point_d& q)
    { assert(p.dimension() == q.dimension());
      dim = p.dimension();
      lower = new NT[dim];
      upper = new NT[dim];
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

    inline int dimension() const { return dim;}
    
    inline NT min_coord(int i) const {
      return lower[i];
    }

    inline NT max_coord(int i) const {
      return upper[i];
    }

  }; // end of class

} // namespace CGAL
#endif // CGAL_ISO_RECTANGLE_D_H

