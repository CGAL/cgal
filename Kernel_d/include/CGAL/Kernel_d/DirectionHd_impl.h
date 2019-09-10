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

#ifndef CGAL_DIRECTIONHD_C
#define CGAL_DIRECTIONHD_C
namespace CGAL {

template <class RT, class LA> 
DirectionHd<RT,LA>::DirectionHd(const VectorHd<RT,LA>& v) : Base(v) {}

template <class RT, class LA>
VectorHd<RT,LA>  DirectionHd<RT,LA>::vector() const
{ return VectorHd<RT,LA>(*this); }

template <class RT, class LA>
DirectionHd<RT,LA>  DirectionHd<RT,LA>::
transform(const Aff_transformationHd<RT,LA>& t) const
{ return vector().transform(t).direction(); }

template <class RT, class LA>
Comparison_result DirectionHd<RT,LA>::
cmp(const DirectionHd<RT,LA>& h1, 
    const DirectionHd<RT,LA>& h2) 
{ 
  if (h1.identical(h2)) return EQUAL; 
  int i, d = h1.dimension(); 
  for (i = 0; i < d && h1.delta(i)==0 && 
              h2.delta(i)==0; i++) {}
  int c1 = CGAL_NTS sign(h1.delta(i)); 
  int c2 = CGAL_NTS sign(h2.delta(i)); 
  if (c1 != c2) return CGAL_NTS compare(c1,c2); 
 
  RT s1 = (RT) CGAL_NTS sign(h2.delta(i)) * h2.delta(i); 
  RT s2 = (RT) CGAL_NTS sign(h1.delta(i)) * h1.delta(i); 

  i++; 
  Comparison_result c; 
  while (i < d) { 
    c = CGAL_NTS compare(s1 * h1.delta(i), s2 * h2.delta(i)); 
    if (c != EQUAL) return c;
    i++; 
  }
  return EQUAL; 
}

template <class RT, class LA>
std::istream& operator>>(std::istream& I, DirectionHd<RT,LA>& dir)
{ dir.copy_on_write(); dir.ptr()->read(I); 
  CGAL_assertion_msg((dir.D()>=0), 
  "operator>>: denominator of direction must be nonnegative."); 
  return I; 
} 

template <class RT, class LA>
std::ostream& operator<<(std::ostream& O, const DirectionHd<RT,LA>& dir)
{ dir.ptr()->print(O,"DirectionHd"); return O; } 


//----------------------- end of file ----------------------------------


} //namespace CGAL
#endif // CGAL_DIRECTIONHD_C
