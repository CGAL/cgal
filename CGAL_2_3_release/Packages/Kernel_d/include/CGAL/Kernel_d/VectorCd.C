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
// file          : include/CGAL/Kernel_d/VectorCd.C
// package       : Kernel_d (0.9.19)
// maintainer    : Michael Seel <seel@mpi-sb.mpg.de>
// author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
// coordinator   : Susan Hert <hert@mpi-sb.mpg.de>
//
// ======================================================================

#ifndef CGAL_VECTORCD_C
#define CGAL_VECTORCD_C
CGAL_BEGIN_NAMESPACE
#define PointCd PointCd2

template <class FT,class LA>
PointCd<FT,LA> VectorCd<FT,LA>::to_point() const
{ return PointCd<FT,LA>(Base(*this)); }

template <class FT,class LA>
PointCd<FT,LA> 
operator+ (const Origin&, const VectorCd<FT,LA>& v)
{ return v.to_point(); }

template <class FT, class LA>
DirectionCd<FT,LA>  VectorCd<FT,LA>::
direction() const
{ CGAL_assertion_msg(!is_zero(), "VectorCd::direction: \
  zero vector cannot be a direction."); 
  return DirectionCd<FT,LA>(*this);
}

template <class FT, class LA>
VectorCd<FT,LA> VectorCd<FT,LA>::
transform(const Aff_transformationCd<FT,LA>& t) const
{ typename LA::Vector res = t.transform_linearly(vector_rep());
  return VectorCd<FT,LA>(dimension(),res.begin(),res.end());
}

template <class FT, class LA>
VectorCd<FT,LA> operator*(const int& n, const VectorCd<FT,LA>& v) 
{ return v.scale(n); }

template <class FT, class LA>
VectorCd<FT,LA> operator*(const FT& n, const VectorCd<FT,LA>& v) 
{ return v.scale(n); }

template <class FT, class LA>
std::istream& operator>>(std::istream& I, VectorCd<FT,LA>& v)
{ v.copy_on_write(); v.ptr->read(I);
  return I; 
}

template <class FT, class LA>
std::ostream& operator<<(std::ostream& O, const VectorCd<FT,LA>& v)
{ v.ptr->print(O,"VectorCd"); return O; } 

template <class FT, class LA>
inline CGAL::io_Operator io_tag(const VectorCd<FT,LA>&) 
{ return CGAL::io_Operator(); }

#undef PointCd
CGAL_END_NAMESPACE
#endif // CGAL_VECTORCD_C

