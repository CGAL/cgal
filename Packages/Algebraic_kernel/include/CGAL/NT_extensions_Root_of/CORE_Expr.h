// Copyright (c) 2003  INRIA Sophia-Antipolis (France) and
//                     Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// Authors : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//           Sylvain Pion     <Sylvain.Pion@sophia.inria.fr>
//           Athanasios Kakargias <grad0460@di.uoa.gr>
//
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473
// (CGAL - Effective Computational Geometry for Curves and Surfaces)

// file : include/CGAL/Root_of/CORE_Expr.h

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
