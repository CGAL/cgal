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
// release       :
// release_date  :
//
// file          : include/CGAL/Plane_separator.h
// package       : ASPAS
// revision      : 3.0
// revision_date : 2003/07/10 
// authors       : Hans Tangelder (<hanst@cs.uu.nl>)
// maintainer    : Hans Tangelder (<hanst@cs.uu.nl>)
// coordinator   : Utrecht University
//
// ======================================================================
#ifndef CGAL_PLANE_SEPARATOR_H
#define CGAL_PLANE_SEPARATOR_H
namespace CGAL {
template < class NT> class Plane_separator {

  public:

  int cutting_dim;
  NT cutting_val;

  inline int cutting_dimension() const { return cutting_dim;}
  inline NT cutting_value() const { return cutting_val;}

  void set_cutting_dimension(int d) {
	  cutting_dim=d;
  }

  void set_cutting_value(NT val) {
	  cutting_val=val;
  }  

  template <class Point>  
  inline bool has_on_negative_side(const Point& i) {
    return i[cutting_dimension()] < cutting_value();
  }

  Plane_separator(const int d, const NT& v) : 
		cutting_dim(d), cutting_val(v) {}
  Plane_separator(const Plane_separator<NT>& s) : 
		cutting_dim(s.cutting_dimension()), 
		cutting_val(s.cutting_value()) {}
  explicit Plane_separator() : cutting_dim(0), cutting_val(0) {}
  Plane_separator<NT>& operator= (const Plane_separator<NT>& s) {
    cutting_dim = s.cutting_dimension();
    cutting_val = s.cutting_value();
    return *this;
  }

  ~Plane_separator() {}

};

 template < class NT> 
 std::ostream& operator<< (std::ostream& s, Plane_separator<NT>& x) {
   s << "\n Separator coordinate: " << x.cutting_dimension() <<
     " value: " << x.cutting_value() << "\n";
   return s;
 }
} // namespace CGAL
#endif // CGAL_PLANE_SEPARATOR_H

