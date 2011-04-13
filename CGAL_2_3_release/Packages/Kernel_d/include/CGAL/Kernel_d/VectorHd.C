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
// file          : include/CGAL/Kernel_d/VectorHd.C
// package       : Kernel_d (0.9.19)
// maintainer    : Michael Seel <seel@mpi-sb.mpg.de>
// author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
// coordinator   : Susan Hert <hert@mpi-sb.mpg.de>
//
// ======================================================================

#ifndef CGAL_VECTORHD_C
#define CGAL_VECTORHD_C
CGAL_BEGIN_NAMESPACE

template <class RT,class LA>
PointHd<RT,LA> VectorHd<RT,LA>::to_point() const
{ return PointHd<RT,LA>(Base(*this)); }

template <class RT,class LA>
PointHd<RT,LA> 
operator+ (const Origin&, const VectorHd<RT,LA>& v)
{ return v.to_point(); }

template <class RT, class LA>
DirectionHd<RT,LA>  VectorHd<RT,LA>::
direction() const
{ CGAL_assertion_msg(!is_zero(), "VectorHd::direction: \
  zero vector cannot be a direction."); 
  return DirectionHd<RT,LA>(*this); 
}

template <class RT, class LA>
VectorHd<RT,LA> VectorHd<RT,LA>::
transform(const Aff_transformationHd<RT,LA>& t) const
{ typename LA::Matrix m_at = t.matrix(); 
  int d = t.dimension(); 
  for (int i = 0; i < d; i++) m_at(i,d) = 0;
  typename LA::Vector res(m_at*vector_rep());
  return VectorHd<RT,LA>(dimension(),res.begin(),res.end()); 
}

template <class RT, class LA>
std::istream& operator>>(std::istream& I, VectorHd<RT,LA>& v)
{ v.copy_on_write(); v.ptr->read(I); 
  CGAL_assertion_msg((v.homogeneous(v.dimension()) > 0),
  "operator>>: denominator of vector must be larger than zero.");
  return I; 
}

template <class RT, class LA>
std::ostream& operator<<(std::ostream& O, const VectorHd<RT,LA>& v)
{ v.ptr->print(O,"VectorHd"); return O; } 

template <class RT, class LA>
inline CGAL::io_Operator io_tag(const VectorHd<RT,LA>&) 
{ return CGAL::io_Operator(); }


CGAL_END_NAMESPACE
#endif // CGAL_VECTORHD_C

