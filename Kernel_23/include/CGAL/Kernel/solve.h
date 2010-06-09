// Copyright (c) 1999  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
  FT denom = b2*c1*a3-b1*c2*a3+c3*b1*a2+b3*c2*a1-c1*b3*a2-b2*c3*a1;

  x = - (b2*c3*d1-b2*c1*d3+c1*b3*d2+b1*c2*d3-c3*b1*d2-b3*c2*d1)/denom;

  z = (b2*d1*a3-b2*a1*d3+b1*a2*d3-b1*d2*a3-d1*b3*a2+a1*b3*d2)/denom;

  y = (a2*c3*d1-a2*c1*d3-c2*d1*a3+c2*a1*d3+d2*c1*a3-d2*c3*a1)/denom;
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
