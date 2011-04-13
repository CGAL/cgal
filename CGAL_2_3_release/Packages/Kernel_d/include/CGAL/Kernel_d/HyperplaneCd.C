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
// file          : include/CGAL/Kernel_d/HyperplaneCd.C
// package       : Kernel_d (0.9.19)
// maintainer    : Michael Seel <seel@mpi-sb.mpg.de>
// author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
// coordinator   : Susan Hert <hert@mpi-sb.mpg.de>
//
// ======================================================================

#ifndef CGAL_HYPERPLANECD_C
#define CGAL_HYPERPLANECD_C
CGAL_BEGIN_NAMESPACE

template <class FT, class LA>
VectorCd<FT,LA>  HyperplaneCd<FT,LA>::
orthogonal_vector() const
{ return VectorCd<FT,LA>(dimension(),coefficients_begin(),
			             coefficients_end()-1); }

template <class FT, class LA>
Comparison_result HyperplaneCd<FT,LA>::
weak_cmp(const HyperplaneCd<FT,LA>& h1,
         const HyperplaneCd<FT,LA>& h2)
{ 
  CGAL_assertion_msg((h1.dimension()==h2.dimension()), 
    "HyperplaneCd::weak_cmp: dimensions disagree.");
  if(h1.identical(h2)) return EQUAL;

  int i, d = h1.dimension();
  for (i = 0; i <= d && 
              h1.coefficient(i) == FT(0) && 
              h2.coefficient(i) == FT(0); i++); // no body
  if (h1.coefficient(i) == FT(0)) return SMALLER;
  if (h2.coefficient(i) == FT(0)) return LARGER;
 
  int s = CGAL_NTS sign(h1.coefficient(i)) * 
          CGAL_NTS sign(h2.coefficient(i));
  FT s1 = (FT)s * h2.coefficient(i);
  FT s2 = (FT)s * h1.coefficient(i);
  // |s1 * h1.coefficient(i)| is 
  // $\Labs{ |h1.coefficient(i)*h2.coefficient(i)| }$

  Comparison_result c;
  while (++i <= d) { 
    c = CGAL_NTS compare(s1 * h1.coefficient(i),
                         s2 * h2.coefficient(i));
    if (c != EQUAL) return c;
  }
  return EQUAL;
}

template <class FT, class LA>
Comparison_result HyperplaneCd<FT,LA>::
strong_cmp(const HyperplaneCd<FT,LA>& h1, 
           const HyperplaneCd<FT,LA>& h2)
{ 
  CGAL_assertion_msg((h1.dimension()==h2.dimension()), 
  "HyperplaneCd::strong_cmp: dimensions disagree.");
  if (h1.identical(h2))  return EQUAL;

  int i;
  int d = h1.dimension();
  for (i = 0; i <=d && 
              h1.coefficient(i)==FT(0) && 
              h2.coefficient(i)==FT(0); i++) ; // no body
  int c1 = CGAL_NTS sign(h1.coefficient(i));
  int c2 = CGAL_NTS sign(h2.coefficient(i));
  if (c1 != c2) return CGAL_NTS compare(c1,c2);
  FT s1 = (FT)CGAL_NTS sign(h2.coefficient(i)) * h2.coefficient(i); 
  FT s2 = (FT)CGAL_NTS sign(h1.coefficient(i)) * h1.coefficient(i);

  Comparison_result c;
  while (++i <= d) { 
    c = CGAL_NTS compare(s1 * h1.coefficient(i), 
                         s2 * h2.coefficient(i));
    if (c != EQUAL) return c;
  }
  return EQUAL;
}

template <class FT, class LA>
std::istream& operator>>(std::istream& I, HyperplaneCd<FT,LA>& h) 
{ h.copy_on_write(); h.ptr->read(I); return I; }

template <class FT, class LA>
std::ostream& operator<<(std::ostream& O, const HyperplaneCd<FT,LA>& h)
{ h.ptr->print(O,"HyperplaneCd"); return O; } 

template <class FT, class LA>
inline CGAL::io_Operator io_tag(const HyperplaneCd<FT,LA>&) 
{ return CGAL::io_Operator(); }
 

CGAL_END_NAMESPACE
#endif // CGAL_HYPERPLANECD_C


