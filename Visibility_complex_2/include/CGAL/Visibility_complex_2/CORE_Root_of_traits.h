// Copyright (c) 2001-2004  ENS of Paris (France).
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
// Author(s)     : Luc Habert


#ifndef CGAL_VISIBILITY_COMPLEX_2_CORE_ROOT_OF_TRAITS_H
#define CGAL_VISIBILITY_COMPLEX_2_CORE_ROOT_OF_TRAITS_H

#include <CGAL/basic.h>
#include<CGAL/CORE_Expr.h>

CGAL_BEGIN_NAMESPACE
namespace Visibility_complex_2_details {

template<class RT> class CORE_Root_of_traits {
  static CORE::Expr convert(const RT&x) {
    return CORE::Expr(CGAL_NTS to_double(x));
  }
public:
  class Root_of_4 {
  public:
    CORE::Expr root;
    Root_of_4 () {};
    Root_of_4 (const RT& c4, const RT& c3, const RT& c2, const RT& c1, 
               const RT& c0, int idx) {
      CORE::Expr ploum[5];
      ploum[0]=convert(c0);
      ploum[1]=convert(c1);
      ploum[2]=convert(c2);
      ploum[3]=convert(c3);
      ploum[4]=convert(c4);
      root=CORE::Expr(CORE::Polynomial<CORE::Expr>(4,ploum),idx+1);
    };
    double to_double() const {
      return root.doubleValue();
    }
  };
  class Root_of_3 {
  public:
    CORE::Expr root;
    Root_of_3 () {};
    Root_of_3 (const RT& c3, const RT& c2, const RT& c1, const RT& c0,
               int idx) {
      CORE::Expr ploum[4];
      ploum[0]=convert(c0);
      ploum[1]=convert(c1);
      ploum[2]=convert(c2);
      ploum[3]=convert(c3);
      root=CORE::Expr(CORE::Polynomial<CORE::Expr>(3,ploum),idx+1);
    };
  };
  class Root_of_2 {
  public:
    CORE::Expr root;
    Root_of_2 () {};
    Root_of_2 (const RT& c2, const RT& c1, const RT& c0, int idx) {
      CORE::Expr ploum[3];
      ploum[0]=convert(c0);
      ploum[1]=convert(c1);
      ploum[2]=convert(c2);
      root=CORE::Expr(CORE::Polynomial<CORE::Expr>(2,ploum),idx+1);
    };
  };
  class Root_of_1 {
 public:
    CORE::Expr root;
    Root_of_1 () {};
    Root_of_1 (const RT& c1, const RT& c0) {
      CORE::Expr ploum[2];
      ploum[0]=convert(c0);
      ploum[1]=convert(c1);
      root=CORE::Expr(CORE::Polynomial<CORE::Expr>(1,ploum),1);
    }
  };
  struct compare_object {
    template<class RO1,class RO2> CGAL::Comparison_result operator () 
      (const RO1& a,const RO2& b) {
      int ploum= a.root.cmp(b.root);
      if (ploum>0) return CGAL::LARGER; else 
        if (ploum==0) return CGAL::EQUAL; else return CGAL::SMALLER;
    };
  };
  class Polynom_1 {
    CORE::Expr coeffs[2];
  public:
    Polynom_1() {}
    Polynom_1(const RT& c1,const RT& c0)  {
      coeffs[0]=convert(c0);
      coeffs[1]=convert(c1);
    }
    template<class RO> CGAL::Sign sign_at(RO &x) {
      CORE::Expr y=coeffs[1]*x.root+coeffs[0];
      int sgn=y.sign();
      if (sgn>0) 
        return POSITIVE; 
      else if (sgn==0)
        return ZERO;
      else
        return NEGATIVE;
    }
  };
  class Polynom_2 {
    CORE::Expr coeffs[3];
  public:
    Polynom_2() {};
    Polynom_2(const RT& c2,const RT& c1,const RT& c0)  {
      coeffs[0]=convert(c0);
      coeffs[1]=convert(c1);
      coeffs[2]=convert(c2);
    }
    template<class RO> CGAL::Sign sign_at(RO& x) {
      CORE::Expr y=(coeffs[2]*x.root+coeffs[1])*x.root+coeffs[0];
      int sgn=y.sign();
      if (sgn>0) 
        return POSITIVE; 
      else if (sgn==0)
        return ZERO;
      else
        return NEGATIVE;
    }
  };
};

}
CGAL_END_NAMESPACE
#endif //CGAL_VISIBILITY_COMPLEX_2_CORE_ROOT_OF_TRAITS_H
