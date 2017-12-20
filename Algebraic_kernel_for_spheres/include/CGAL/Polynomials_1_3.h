// Copyright (c) 2005-2006  INRIA Sophia-Antipolis (France).
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
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s) : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//             Sylvain Pion
//             Pedro Machado
//             Julien Hazebrouck
//             Damien Leroy

#ifndef CGAL_ALGEBRAIC_KERNEL_POLYNOMIALS_1_3_H
#define CGAL_ALGEBRAIC_KERNEL_POLYNOMIALS_1_3_H

#include <CGAL/license/Circular_kernel_3.h>


#include <CGAL/enum.h>

namespace CGAL {

// polynomials of the form aX + +bY + cZ + d
template < typename FT_ >
class Polynomial_1_3
{
  FT_ rep[4]; // stores a, b, c, d
  
public:
  
  typedef FT_ FT;
  
  Polynomial_1_3(){}
  
  Polynomial_1_3(const FT & a, const FT & b, const FT & c, const FT & d)
  { 
    rep[0]=a;
    rep[1]=b;
    rep[2]=c;
    rep[3]=d;
  }

  const FT & a() const
  { return rep[0]; }

  const FT & b() const
  { return rep[1]; }
  
  const FT & c() const
  { return rep[2]; }
  
  const FT & d() const
  { return rep[3]; }

  bool undefined() const {
    return is_zero(a()) &&
           is_zero(b()) &&
           is_zero(c()) &&
           is_zero(d());
  }

  bool empty_space() const {
    return is_zero(a()) &&
           is_zero(b()) &&
           is_zero(c()) &&
           (!is_zero(d()));
  }
};

template < typename FT >
inline
bool 
operator == ( const Polynomial_1_3<FT> & p1,
	      const Polynomial_1_3<FT> & p2 )
{
  return( (p1.a() == p2.a()) && 
	  (p1.b() == p2.b()) &&
	  (p1.c() == p2.c()) &&
	  (p1.d() == p2.d()) );
}

template < typename FT >
inline
bool 
same_solutions ( const Polynomial_1_3<FT> & p1,
	         const Polynomial_1_3<FT> & p2 )
{
  if(p1.undefined()) return p2.undefined();
  if(p1.empty_space()) return p2.empty_space();
  if(p2.undefined()) return false;
  if(p2.empty_space()) return false;
  if(is_zero(p1.a())) {
    if(!is_zero(p2.a())) return false;
    if(is_zero(p1.b())) {
      if(!is_zero(p2.b())) return false;
      return p1.c() * p2.d() == p1.d() * p2.c();
    }
    return (p2.c() * p1.b() == p1.c() * p2.b()) &&
           (p2.d() * p1.b() == p1.d() * p2.b());
  }
  return (p2.b() * p1.a() == p1.b() * p2.a()) &&
         (p2.c() * p1.a() == p1.c() * p2.a()) &&
         (p2.d() * p1.a() == p1.d() * p2.a());
}

} //namespace CGAL

#endif //CGAL_ALGEBRAIC_KERNEL_POLYNOMIALS_1_3_H
