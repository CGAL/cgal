// ======================================================================
//
// Copyright (c) 2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.3-I-67 $
// release_date  : $CGAL_Date: 2001/05/29 $
//
// file          : include/CGAL/Kernel_d/PointHd.C
// package       : Kernel_d (0.9.19)
// maintainer    : Michael Seel <seel@mpi-sb.mpg.de>
// author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
// coordinator   : Susan Hert <hert@mpi-sb.mpg.de>
//
// ======================================================================

#ifndef CGAL_POINTHD_C
#define CGAL_POINTHD_C
CGAL_BEGIN_NAMESPACE

template <class RT, class LA>
PointHd<RT,LA> PointHd<RT,LA>::
transform(const Aff_transformationHd<RT,LA>& t) const
{ typename LA::Vector res = t(vector_rep());
  return PointHd<RT,LA>(dimension(),res.begin(),res.end()); }

template <class RT, class LA>
VectorHd<RT,LA> PointHd<RT,LA>::operator-(const Origin&) const 
{ return VectorHd<RT,LA>(Base(*this)); }

template <class RT, class LA>
PointHd<RT,LA> PointHd<RT,LA>::operator+(const VectorHd<RT,LA> &v) const
{ PointHd<RT,LA> res(dimension()); 
  res.ptr->homogeneous_add(ptr, v.ptr); 
  return res; 
}

template <class RT, class LA>
PointHd<RT,LA> PointHd<RT,LA>::operator-(const VectorHd<RT,LA> &v) const
{ PointHd<RT,LA> res(dimension()); 
  res.ptr->homogeneous_sub(ptr, v.ptr); 
  return res; 
}

template <class RT, class LA>
PointHd<RT,LA>& PointHd<RT,LA>::operator+= (const VectorHd<RT,LA>& v)
{ int d = dimension(); 
  PointHd<RT,LA> old(*this); 
  *this = PointHd<RT,LA>(d); 
  ptr->homogeneous_add(old.ptr, v.ptr); 
  return *this; 
}

template <class RT, class LA>
PointHd<RT,LA>& PointHd<RT,LA>::operator-= (const VectorHd<RT,LA>& v)
{ int d = dimension(); 
  PointHd<RT,LA> old(*this); 
  *this = PointHd<RT,LA>(d); 
  ptr->homogeneous_sub(old.ptr, v.ptr); 
  return *this; 
}

template <class RT, class LA>
std::istream& operator>>(std::istream& I, PointHd<RT,LA>& p)
{ p.copy_on_write(); p.ptr->read(I); 
  CGAL_assertion_msg((p.homogeneous(p.dimension()) > 0), 
  "operator>>: denominator of point must be larger than zero.");
  return I; 
}

template <class RT, class LA>
std::ostream& operator<<(std::ostream& O, const PointHd<RT,LA>& p)
{ p.ptr->print(O,"PointHd"); return O; } 

template <class RT, class LA>
inline CGAL::io_Operator io_tag(const PointHd<RT,LA>&) 
{ return CGAL::io_Operator(); }



CGAL_END_NAMESPACE
#endif // CGAL_POINTHD_C

