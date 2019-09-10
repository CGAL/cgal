// Copyright (c) 2001  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Michael Seel <seel@mpi-sb.mpg.de>

#ifndef CGAL_HYPERPLANECD_C
#define CGAL_HYPERPLANECD_C
namespace CGAL {

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
              h2.coefficient(i) == FT(0); i++) {}
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
              h2.coefficient(i)==FT(0); i++) {}
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
{ h.copy_on_write(); h.ptr()->read(I); return I; }

template <class FT, class LA>
std::ostream& operator<<(std::ostream& O, const HyperplaneCd<FT,LA>& h)
{ h.ptr()->print(O,"HyperplaneCd"); return O; } 

} //namespace CGAL
#endif // CGAL_HYPERPLANECD_C
