// Copyright (c) 2003-2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Monique Teillaud, Sylvain Pion

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_ALGEBRAIC_KERNEL_FOR_CIRCLES_POLYNOMIALS_2_2_H
#define CGAL_ALGEBRAIC_KERNEL_FOR_CIRCLES_POLYNOMIALS_2_2_H

//////////// FIXME - pb RT (cas general Polynomial_2_2) ou FT (ici)

#include <CGAL/enum.h>

namespace CGAL {

// polynomials of the form (X-a)^2 + (Y-b)^2 - R^2
template < typename FT_ >
class Polynomial_for_circles_2_2
{
  FT_ rep[3]; // stores a, b, R^2
  
public:
  
  typedef FT_ FT;
  
  Polynomial_for_circles_2_2(){}
  
  Polynomial_for_circles_2_2(const FT & a, const FT & b, const FT & rsq)
  { 
    rep[0]=a;
    rep[1]=b;
    rep[2]=rsq;
  }

  const FT & a() const
  { return rep[0]; }

  const FT & b() const
  { return rep[1]; }
  
  const FT & r_sq() const
  { return rep[2]; }
};

template < typename FT >
bool 
operator == ( const Polynomial_for_circles_2_2<FT> & p1,
	      const Polynomial_for_circles_2_2<FT> & p2 )
{
  return( (p1.a() == p2.a()) && 
              (p1.b() == p2.b()) &&
              (p1.r_sq() == p2.r_sq()) );
}

} //namespace CGAL

#endif //CGAL_ALGEBRAIC_KERNEL_FOR_CIRCLES_POLYNOMIALS_2_2_H
