// Copyright (c) 1999  
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
// 
//
// Author(s)     : Andreas Fabri
 

#ifndef CGAL_SOLVE_H
#define CGAL_SOLVE_H

namespace CGAL {


template <class FT>
void solve (const FT &a1, const FT &a2, const FT &a3,
            const FT &b1, const FT &b2, const FT &b3,
            const FT &c1, const FT &c2, const FT &c3,
            const FT &d1, const FT &d2, const FT &d3,
            FT &x, FT &y, FT &z)
{
#if 0
  FT denom = CGAL::determinant(a1, a2, a3, b1, b2, b3, c1, c2, c3);

  x = CGAL::determinant(b1, b2, b3, c1, c2, c3, d1, d2, d3)/denom;

  y = - CGAL::determinant(a1, a2, a3, c1, c2, c3, d1, d2, d3)/denom;

  z = CGAL::determinant(a1, a2, a3, b1, b2, b3, d1, d2, d3)/denom;
#else
  // Same as above, but expanded to factorize 6 internal 2x2 minors.

  FT ab23 = a3*b2 - a2*b3;
  FT ab13 = a3*b1 - a1*b3;
  FT ab12 = a2*b1 - a1*b2;

  FT denom = ab23*c1 - ab13*c2 + ab12*c3;

  FT cd23 = c3*d2 - c2*d3;
  FT cd13 = c3*d1 - c1*d3;
  FT cd12 = c2*d1 - c1*d2;

  x = (b3*cd12 - b2*cd13 + b1*cd23)/denom;

  y = (a2*cd13 - cd12*a3 - cd23*a1)/denom;

  z = (ab23*d1 + ab12*d3 - ab13*d2)/denom;
#endif
}


// this is for a parabola c1, c2, c3 are equal to 1
template <class FT>
void solve_quadratic (const FT &a1, const FT &a2, const FT &a3,
                      const FT &b1, const FT &b2, const FT &b3,
                      const FT &d1, const FT &d2, const FT &d3,
                      FT &x, FT &y, FT &z)
{
  FT denom = b2*a3-b1*a3+b1*a2+b3*a1-b3*a2-b2*a1;

  x = - (b2*d1-b2*d3+b3*d2+b1*d3-b1*d2-b3*d1)/denom;

  z = (b2*d1*a3-b2*a1*d3+b1*a2*d3-b1*d2*a3-d1*b3*a2+a1*b3*d2)/denom;

  y = (a2*d1-a2*d3-d1*a3+a1*d3+d2*a3-d2*a1)/denom;
}


} //namespace CGAL

#endif // CGAL_SOLVE_H
