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
// file          : include/CGAL/Kernel_d/DirectionCd.C
// package       : Kernel_d (0.9.19)
// maintainer    : Michael Seel <seel@mpi-sb.mpg.de>
// author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
// coordinator   : Susan Hert <hert@mpi-sb.mpg.de>
//
// ======================================================================

#ifndef CGAL_DIRECTIONCD_C
#define CGAL_DIRECTIONCD_C
CGAL_BEGIN_NAMESPACE

template <class FT, class LA> 
DirectionCd<FT,LA>::DirectionCd(const VectorCd<FT,LA>& v) : Base(v) {}

template <class FT, class LA>
VectorCd<FT,LA> DirectionCd<FT,LA>::vector() const
{ return VectorCd<FT,LA>(*this); }

template <class FT, class LA>
DirectionCd<FT,LA>  DirectionCd<FT,LA>::
transform(const Aff_transformationCd<FT,LA>& t) const
{ return vector().transform(t).direction(); }

template <class FT, class LA>
Comparison_result DirectionCd<FT,LA>::
cmp(const DirectionCd<FT,LA>& h1, 
    const DirectionCd<FT,LA>& h2) 
{ 
  if (h1.identical(h2)) return EQUAL; 
  int i, d = h1.dimension(); 
  for (i = 0; i < d && h1.delta(i)==FT(0) && 
                       h2.delta(i)==FT(0); ++i) ; // no body
  int c1 = CGAL_NTS sign(h1.delta(i));
  int c2 = CGAL_NTS sign(h2.delta(i));
  if (c1 != c2) return CGAL_NTS compare(c1,c2);
 
  FT s1 = (FT) CGAL_NTS sign(h2.delta(i)) * h2.delta(i); 
  FT s2 = (FT) CGAL_NTS sign(h1.delta(i)) * h1.delta(i); 

  i++;
  Comparison_result c; 
  while (i < d) { 
    c = CGAL_NTS compare(s1 * h1.delta(i), s2 * h2.delta(i));
    if (c != EQUAL) return c;
    i++;
  }
  return EQUAL;
}

template <class FT, class LA>
std::istream& operator>>(std::istream& I, DirectionCd<FT,LA>& dir)
{ dir.copy_on_write(); dir.ptr->read(I);
  return I; 
} 

template <class FT, class LA>
std::ostream& operator<<(std::ostream& O, const DirectionCd<FT,LA>& dir)
{ dir.ptr->print(O,"DirectionCd"); return O; } 

template <class FT, class LA>
inline CGAL::io_Operator io_tag(const DirectionCd<FT,LA>&) 
{ return CGAL::io_Operator(); }

CGAL_END_NAMESPACE
#endif // CGAL_DIRECTIONCD_C
//----------------------- end of file ----------------------------------


