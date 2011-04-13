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
// file          : include/CGAL/Kernel_d/PointCd.C
// package       : Kernel_d (0.9.19)
// maintainer    : Michael Seel <seel@mpi-sb.mpg.de>
// author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
// coordinator   : Susan Hert <hert@mpi-sb.mpg.de>
//
// ======================================================================

#ifndef CGAL_POINTCDXXX_C
#define CGAL_POINTCDXXX_C
CGAL_BEGIN_NAMESPACE
#define PointCd PointCd2

template <class FT, class LA>
PointCd<FT,LA> PointCd<FT,LA>::
transform(const Aff_transformationCd<FT,LA>& t) const
{ typename LA::Vector res = t(vector_rep());
  return PointCd<FT,LA>(dimension(),res.begin(),res.end()); }

template <class FT, class LA>
VectorCd<FT,LA> PointCd<FT,LA>::operator-(const Origin&) const 
{ return VectorCd<FT,LA>(Base(*this)); }

template <class FT, class LA>
PointCd<FT,LA> PointCd<FT,LA>::operator+(const VectorCd<FT,LA> &v) const
{ PointCd<FT,LA> res(dimension()); 
  res.ptr->cartesian_add(ptr, v.ptr);
  return res; 
}

template <class FT, class LA>
PointCd<FT,LA> PointCd<FT,LA>::operator-(const VectorCd<FT,LA> &v) const
{ PointCd<FT,LA> res(dimension()); 
  res.ptr->cartesian_sub(ptr, v.ptr); 
  return res; 
}

template <class FT, class LA>
PointCd<FT,LA>& PointCd<FT,LA>::operator+= (const VectorCd<FT,LA>& v)
{ int d = dimension(); 
  PointCd<FT,LA> old(*this); 
  *this = PointCd<FT,LA>(d); 
  ptr->cartesian_add(old.ptr, v.ptr); 
  return *this; 
}

template <class FT, class LA>
PointCd<FT,LA>& PointCd<FT,LA>::operator-= (const VectorCd<FT,LA>& v)
{ int d = dimension(); 
  PointCd<FT,LA> old(*this); 
  *this = PointCd<FT,LA>(d); 
  ptr->cartesian_sub(old.ptr, v.ptr); 
  return *this; 
}

template <class FT, class LA>
std::istream& operator>>(std::istream& I, PointCd<FT,LA>& p)
{ p.copy_on_write(); p.ptr->read(I);
  return I; 
}

template <class FT, class LA>
std::ostream& operator<<(std::ostream& O, const PointCd<FT,LA>& p)
{ p.ptr->print(O,"PointCd"); return O; } 

template <class FT, class LA>
inline CGAL::io_Operator io_tag(const PointCd<FT,LA>&) 
{ return CGAL::io_Operator(); }

#undef PointCd
CGAL_END_NAMESPACE
#endif // CGAL_POINTCDXXX_C

