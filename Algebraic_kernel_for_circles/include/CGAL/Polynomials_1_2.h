// Copyright (c) 2003-2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Monique Teillaud, Sylvain Pion

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473
// (ECG - Effective Computational Geometry for Curves and Surfaces)
// and a STREP (FET Open) Project under Contract No  IST-006413
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_ALGEBRAIC_KERNEL_FOR_CIRCLES_POLYNOMIALS_1_2_H
#define CGAL_ALGEBRAIC_KERNEL_FOR_CIRCLES_POLYNOMIALS_1_2_H

#include <CGAL/license/Circular_kernel_2.h>


#include <CGAL/enum.h>

namespace CGAL {

template < typename RT_ >
class Polynomial_1_2
{
  RT_ rep[3]; // stores a, b, c for line ax+by+c=0

public:

  typedef RT_ RT;

  Polynomial_1_2(){}

  Polynomial_1_2(const RT & a, const RT & b, const RT & c)
  {
    rep[0]=a;
    rep[1]=b;
    rep[2]=c;
  }

  const RT & a() const
  { return rep[0]; }

  const RT & b() const
  { return rep[1]; }

  const RT & c() const
  { return rep[2]; }
};

template < typename RT >
bool
operator == ( const Polynomial_1_2<RT> & p1,
              const Polynomial_1_2<RT> & p2 )
{
  return( (p1.a() == p2.a()) &&
              (p1.b() == p2.b()) &&
              (p1.c() == p2.c()) );
}

} //namespace CGAL

#endif //CGAL_ALGEBRAIC_KERNEL_FOR_CIRCLES_POLYNOMIALS_1_2_H
