// Copyright (c) 2003-2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Monique Teillaud, Sylvain Pion, Athanasios Kakargias

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_ROOT_OF_CORE_EXPR_H
#define CGAL_ROOT_OF_CORE_EXPR_H

#include <CGAL/CORE_Expr.h>
#include <CGAL/Root_of/Root_of_traits.h>

#ifdef CORE_USE_GMPXX
#  warning Core uses gmpxx
#  include <CGAL/gmpxx.h>
#endif

namespace CGAL {

  namespace CGALi {

    struct Core_traits
    {
      typedef CORE::Expr RootOf_1;
      typedef CORE::Expr RootOf_2;
      typedef CORE::Expr RootOf_3;
      typedef CORE::Expr RootOf_4;
    };

  } // namespace CGALi

    //Coeffs are Core::BigInt

    inline
    CORE::Expr
    make_root_of_2(const CORE::BigInt &a, const CORE::BigInt &b,
                   const CORE::BigInt &c, bool d)
    {
      CORE::BigInt discriminant = b*b - 4*a*c;
      assert( discriminant >= 0);
      if (discriminant == 0)
      {
        return CORE::Expr(CORE::BigRat(-b,2*a));
      }
      else
      {
        CORE::BigInt coeffs[3] = {c,b,a};
        Polynomial< CORE::BigInt > p(2,coeffs); // degree,coeffsarray[degree+1]
        // 1 for the smaller 2 for the larger
        return rootOf(p, d ? 1 : 2 );
      }
    }

    template <>
    struct Root_of_traits< CORE::BigInt > 
      : public CGALi::Core_traits {};

    //General Case Coeffs are Core::Expr
#ifndef CORE_USE_GMPXX
    inline
    CORE::Expr
    make_root_of_2(const CORE::Expr &a, const CORE::Expr &b,
                   const CORE::Expr &c, bool d)
    {
      return CGALi::make_root_of_2_sqrt(a,b,c,d);
    }

    template <>
    struct Root_of_traits< CORE::Expr >
      : public CGALi::Core_traits {};
#endif

}

#endif // CGAL_ROOT_OF_CORE_EXPR_H
