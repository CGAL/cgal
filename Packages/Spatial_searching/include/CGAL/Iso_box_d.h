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
#include <CGAL/representation_tags.h>
#include <functional>
#include <algorithm>
#include <new>
#include <cassert>

namespace CGAL {
  
  template <class Kernel> class Iso_box_d {
  public:
    typedef typename Kernel::FT FT;
    typedef typename Kernel::Point_d Point_d;
    typedef typename Kernel::Rep_tag Rep_tag;
   
  private:

    int dim;
    Point_d* lower;
    Point_d* upper;

  
  void construct_box(const Point_d& p, const Point_d& q, Cartesian_tag tag) {
     FT low[dim];
     FT upp[dim];
     for (int i = 0; i < dim; ++i) {
	  if (p[i] <= q[i]) {
		low[i]=p[i]; 
                upp[i]=q[i];
	  }
	  else {
		low[i]=q[i]; 
                upp[i]=p[i];
	  }
     }	  
     lower = new Point_d(dim,low,low+dim);
     upper = new Point_d(dim,upp,upp+dim);
}  
  
void construct_box(const Point_d& p, const Point_d& q, Homogeneous_tag tag) {
 {
     typedef typename Kernel::RT RT;
     RT low[dim];
     RT upp[dim];
     for (int i = 0; i < dim; ++i) {
	  if (p[i] <= q[i]) {
		low[i]=p.homogeneous(i); 
                upp[i]=q.homogeneous(i);
	  }
	  else {
		low[i]=q.homogeneous(i); 
                upp[i]=p.homogeneous(i);
	  }
     }	  
     lower = new Point_d(dim,low,low+dim,RT(1));
     upper = new Point_d(dim,upp,upp+dim,RT(1));
  }
}

  public:

    Iso_box_d(const Point_d& p, const Point_d& q)
    { CGAL_precondition(p.dimension() == q.dimension());
      dim = p.dimension();
      typename Kernel::Rep_tag tag;
#if defined(__sun) && defined(__SUNPRO_CC)
     // to avoid a warning "tag has not yet been assigned a value"
     // typedef typename Kernel::Rep_tag Rep_tag;
     tag = Rep_tag();
#endif // SUNPRO
     construct_box(p,q,tag);  
    }
  
    // copy constructor
    Iso_box_d(const Iso_box_d& b) : dim(b.dim) {
      lower = new Point_d(*(b.lower));
      upper = new Point_d(*(b.upper));
    }


    // destructor
    ~Iso_box_d() {
	delete lower;
	delete upper;
    }


    Bounded_side bounded_side(const Point_d& p) const
    { 
      CGAL_precondition(p.dimension() == dim);
      FT h;
      for (int i = 0; i < dim; ++i) {
        h=p[i];
        if ( (h < (*lower)[i]) || (h > (*upper)[i]) ) return ON_UNBOUNDED_SIDE;
      }
      for (int i = 0; i < dim; ++i) {
        h=p[i];
        if ( (h == (*lower)[i]) || (h == (*upper)[i]) ) return ON_BOUNDARY;
      }
      return ON_BOUNDED_SIDE;
    }

    bool has_on_bounded_side(const Point_d& p) const
    { 
      return (bounded_side(p)==ON_BOUNDED_SIDE);
    } 

    bool has_on_unbounded_side(const Point_d& p) const
    {
      return (bounded_side(p)==ON_UNBOUNDED_SIDE); 
    } 

    bool has_on_boundary(const Point_d& p) const
    {
      return (bounded_side(p)==ON_BOUNDARY); 
    } 

      
    int dimension() const { return dim;}
    
    FT min_coord(int i) const {
      return (*lower)[i];
    }

    FT max_coord(int i) const {
      return (*upper)[i];
    }

   Point_d min() const
   {
     return Point_d(*lower);
   }

   Point_d max() const
   {
     return Point_d(*upper);
   }
 
   
  FT volume() const
  {
     FT  vol = (*upper)[0]-(*lower)[0];
     for (int i = 1; i < dim; ++i) {
        vol=vol*((*upper)[i]-(*lower)[i]);
     }
     return vol;
  }

bool is_degenerate() const
  {
     for (int i = 0; i < dim; ++i) {
        if ((*lower)[i]==(*upper)[i]) return true;
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

