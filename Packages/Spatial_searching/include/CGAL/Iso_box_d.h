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
  
  template <class Kernel> class Iso_box_d {
  public:
    typedef typename Kernel::FT FT;
    typedef typename Kernel::Point_d Point_d;
   
  private:

    int dim;
    FT *lower;
    FT *upper;
    
  public:

    Iso_box_d(const Point_d& p, const Point_d& q)
    { CGAL_precondition(p.dimension() == q.dimension());
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
      lower = new FT[dim];
      upper = new FT[dim];
      for (int i = 0; i < dim; ++i) {
		lower[i]=b.lower[i]; 
                upper[i]=b.upper[i];
      }
    }


    // destructor
    ~Iso_box_d() {
	delete [] lower;
	delete [] upper;
    }

    bool has_on_bounded_side(const Point_d& p) const
    { 
      CGAL_precondition(p.dimension() == dim);
      FT h;
      for (int i = 0; i < dim; ++i) {
        h=p[i];
        if ( (h <= lower[i]) || (h >= upper[i]) ) return false;
      }
      return true;
    } 

    bool has_on_unbounded_side(const Point_d& p) const
    {
      CGAL_precondition(p.dimension() == dim);
      FT h;
      for (int i = 0; i < dim; ++i) {
        h=p[i];
        if ( (h >= lower[i]) && (h <= upper[i]) ) return false;
      }
      return true;
    } 

    bool has_on_boundary(const Point_d& p) const
    {
      CGAL_precondition(p.dimension() == dim);
      if (has_on_unbounded_side(p)) return false; 
      FT h;
      for (int i = 0; i < dim; ++i) {
        h=p[i];
        if ( (h == lower[i]) || (h == upper[i]) ) return true;
      }
      return false;
    } 

    Bounded_side bounded_side(const Point_d& p) const
    { 
      CGAL_precondition(p.dimension() == dim);
      if (has_on_unbounded_side(p)) return ON_UNBOUNDED_SIDE; 
      FT h;
      for (int i = 0; i < dim; ++i) {
        h=p[i];
        if ( (h == lower[i]) || (h == upper[i]) ) return ON_BOUNDARY;
      }
      return ON_BOUNDED_SIDE;
    }
      
    int dimension() const { return dim;}
    
    FT min_coord(int i) const {
      return lower[i];
    }

    FT max_coord(int i) const {
      return upper[i];
    }

/* does not work, because assignment to p[i] is not allowed
   Point_d min() const
   {
     Point_d p(dim,ORIGIN);
     for (int i = 0; i < dim; ++i) {
        p[i]=lower[i];
     }
     return Point_d(p);
  }

  Point_d max() const
  {
     Point_d p(dim,ORIGIN);
     for (int i = 0; i < dim; ++i) {
        p[i]=upper[i];
     }
     return Point_d(p);
  }
 */
  /* use of Base_vector() after Kernel_d\interface-test.C 
  Point_d min() const
   {
     Point_d p(dim,ORIGIN);
     for (int i = 0; i < dim; ++i) {
        // v is ith base vector
        Vector_d v = Vector_d(dim,Vector_d::Base_vector(),i);
        p=p+lower[i]*v;
     }
     return Point_d(p);
  }

  Point_d max() const
  {
     Point_d p(dim,ORIGIN);
     for (int i = 0; i < dim; ++i) {
        // v is ith base vector
        Vector_d v = Vector_d(dim,Vector_d::Base_vector(),i);
        p=p+upper[i]*v;
     }
     return Point_d(p);
  }*/

  FT volume() const
  {
     FT  vol = upper[0]-lower[0];
     for (int i = 1; i < dim; ++i) {
        vol=vol*(upper[i]-lower[i]);
     }
     return vol;
  }

bool is_degenerate() const
  {
     for (int i = 0; i < dim; ++i) {
        if (lower[i]==upper[i]) return true;
     }
     return false;
  }

}; // end of class

template <class Kernel>
inline bool
operator!=(const Iso_box_d<Kernel>& b1, Iso_box_d<Kernel>& b2)
{
  CGAL_precondition(b1.dimension() == b2.dimension());
  unsigned int dim = b1.dimension();
  for (unsigned int i = 0; i < dim; ++i) {
			if ( (b1.min_coord(i) != b2.min_coord(i)) || (b1.max_coord(i) != b2.max_coord(i)) ) return true;
		}
  return false; 
}

template <class Kernel>
inline bool
operator==(const Iso_box_d<Kernel>& b1, Iso_box_d<Kernel>& b2)
{
  CGAL_precondition(b1.dimension() == b2.dimension());
  unsigned int dim = b1.dimension();
  for (unsigned int i = 0; i < dim; ++i) {
			if ( (b1.min_coord(i) != b2.min_coord(i)) || (b1.max_coord(i) != b2.max_coord(i)) ) return false;
		}
  return true; 
}
  


} // namespace CGAL
#endif // CGAL_ISO_BOX_D_H

