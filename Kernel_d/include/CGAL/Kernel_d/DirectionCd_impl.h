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

#ifndef CGAL_DIRECTIONCD_C
#define CGAL_DIRECTIONCD_C
namespace CGAL {

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
                       h2.delta(i)==FT(0); ++i) {}
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
{ dir.copy_on_write(); dir.ptr()->read(I);
  return I; 
} 

template <class FT, class LA>
std::ostream& operator<<(std::ostream& O, const DirectionCd<FT,LA>& dir)
{ dir.ptr()->print(O,"DirectionCd"); return O; } 

} //namespace CGAL
#endif // CGAL_DIRECTIONCD_C
//----------------------- end of file ----------------------------------
