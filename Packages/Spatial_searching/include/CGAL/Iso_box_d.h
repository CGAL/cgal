// ======================================================================
//
// Copyright (c) 2002 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.5-I-99 $
// release_date  : $CGAL_Date: 2003/05/23 $
//
// file          : include/CGAL/Iso_box_d.h
// package       : ASPAS (3.12)
// maintainer    : Hans Tangelder <hanst@cs.uu.nl>
// revision      : 2.4 
// revision_date : 2003/02/01 
// authors       : Hans Tangelder (<hanst@cs.uu.nl>)
// coordinator   : Utrecht University
//
// ======================================================================

#ifndef CGAL_ISO_BOX_D_H
#define CGAL_ISO_BOX_D_H
#include <CGAL/enum.h>
#include <functional>
#include <algorithm>
#include <new>
#include <cassert>

namespace CGAL {
  
  template <class R> class Iso_box_d {
  public:
    typedef typename R::FT FT;
    typedef typename R::Point_d Point_d;

  private:

    int dim;
    FT *lower;
    FT *upper;
    
  public:

    Iso_box_d(const Point_d& p, const Point_d& q)
    { assert(p.dimension() == q.dimension());
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
    Iso_box_d(const Iso_box_d& b) : dim(b.dim) {
      std::cout << "dim=" << dim << std::endl;
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
      for (int i = 0; i < dim; ++i) {
        h=p[i];
        if ( (h <= lower[i]) || (h >= upper[i]) ) return false;
      }
      return true;
    } 

    bool has_on_unbounded_side(const Point_d& p) const
    {
      FT h;
      for (int i = 0; i < dim; ++i) {
        h=p[i];
        if ( (h >= lower[i]) && (h <= upper[i]) ) return false;
      }
      return true;
    } 

    bool has_on_boundary(const Point_d& p) const
    {
      return ( !(has_on_unbounded_side(p)) && !(has_on_bounded_side(p)) );
    } 

    Bounded_side bounded_side(const Point_d& p) const
    { 
      if (has_on_bounded_side(p)) return ON_BOUNDED_SIDE;
      if (has_on_unbounded_side(p)) return ON_UNBOUNDED_SIDE;
      return ON_BOUNDARY;
    }
      
    int dimension() const { return dim;}
    
    FT min_coord(int i) const {
      return lower[i];
    }

    FT max_coord(int i) const {
      return upper[i];
    }

    Point_d  operator[](int i) const
    {
    	return  vertex(i);
    }

    Point_d vertex(int i) const
    {
	p = new Point_d(dim,ORIGIN);
	i == i % (2**dim);
        for (int d = dim; d > 0; --d) {
		if (i >= 2**(d-1)) 
			p[d]=upper[i];
		else
			p[d]=lower[i];	
		i = i - 2**(d-1);
        }
	return p;
   }

   Point_d min() const
   {
     p = new Point_d(dim,ORIGIN);
     for (int i = 0; i < dim; ++i) {
        p[i]=lower[i];
     }
     return p;
  }

  Point_d max() const
  {
     p = new Point_d(dim,ORIGIN);
     for (int i = 0; i < dim; ++i) {
        p[i]=upper[i];
     }
     return p;
  }
     
  FT volume() const
  {
     FT  vol = upper[0]-lower[0];
     for (int i = 1; i < dim; ++i) {
        vol=vol*(upper[i]-lower[i]);
     }
     return vol;
  }



inline
bool
operator!=(const Iso_box_d<R>& b)
{
  for (unsigned int i = 0; i < dim; ++i) {
			if ( (lower[i] != b.lower[i]) || (upper[i] != b.upper[i]) ) return true;
		}
  return false; 
}

inline
bool
operator==(const Iso_box_d<R>& b)
{
  for (unsigned int i = 0; i < dim; ++i) {
			if ( (lower[i] != b.lower[i]) || (upper[i] != b.upper[i]) ) return false;
		}
  return true; 
}
  
}; // end of class

} // namespace CGAL
#endif // CGAL_ISO_BOX_D_H

