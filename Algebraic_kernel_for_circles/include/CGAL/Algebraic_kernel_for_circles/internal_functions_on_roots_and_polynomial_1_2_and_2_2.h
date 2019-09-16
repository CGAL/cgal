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
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Monique Teillaud, Sylvain Pion, Julien Hazebrouck

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_ALGEBRAIC_KERNEL_FOR_CIRCLES_FUNCTIONS_ON_ROOTS_AND_POLYNOMIAL_1_2_AND_2_2_H
#define CGAL_ALGEBRAIC_KERNEL_FOR_CIRCLES_FUNCTIONS_ON_ROOTS_AND_POLYNOMIAL_1_2_AND_2_2_H

#include <CGAL/license/Circular_kernel_2.h>



namespace CGAL {
  namespace AlgebraicFunctors {


  template < class AK, class OutputIterator >
  inline 
  OutputIterator
  solve( const typename AK::Polynomial_1_2 & e1,
	 const typename AK::Polynomial_for_circles_2_2 & e2,
	 OutputIterator res )
  {
    typedef typename AK::FT FT;
    typedef typename AK::Root_of_2 Root_of_2;
    typedef typename AK::Root_for_circles_2_2 Root_for_circles_2_2;
    if (is_zero(e1.a())){//horizontal line
      
      const FT hy = -e1.c()/e1.b();
      const FT hdisc = e2.r_sq() - CGAL::square(hy - e2.b());
      CGAL::Sign sign_hdisc = CGAL::sign(hdisc);
    
      if(sign_hdisc == NEGATIVE) return res;
      if(sign_hdisc == ZERO) {	
	*res++ = std::make_pair
	  ( Root_for_circles_2_2(Root_of_2(e2.a()), 
                                 Root_of_2(hy)), 2u);
	return res;
      }
      const Root_of_2 x_res1 = make_root_of_2(e2.a(),FT(-1),hdisc);
      const Root_of_2 x_res2 = make_root_of_2(e2.a(),FT(1),hdisc);
      const Root_of_2 y_res = Root_of_2(hy);  
      *res++ = std::make_pair
	( Root_for_circles_2_2(x_res1, y_res), 1u);
      *res++ = std::make_pair
	( Root_for_circles_2_2(x_res2, y_res), 1u);
      return res;
    }
    else if(is_zero(e1.b())){//vertical line
      
      const FT vx = -e1.c()/e1.a();
      const FT vdisc = e2.r_sq() - CGAL::square(vx - e2.a());
      CGAL::Sign sign_vdisc = CGAL::sign(vdisc);

      if(sign_vdisc == NEGATIVE) return res;
      if(sign_vdisc == ZERO) {
	*res++ = std::make_pair
	  ( Root_for_circles_2_2(Root_of_2(vx), 
                                 Root_of_2(e2.b())), 2u);
	return res;
      }
      
      const Root_of_2 x_res = Root_of_2(vx); 
      const Root_of_2 y_res1 = make_root_of_2(e2.b(),FT(-1),vdisc);
      const Root_of_2 y_res2 = make_root_of_2(e2.b(),FT(1),vdisc);

      *res++ = std::make_pair
	( Root_for_circles_2_2(x_res, y_res1), 1u);
      *res++ = std::make_pair
	( Root_for_circles_2_2(x_res, y_res2), 1u);
      return res;
    }
    else {
      
      const FT line_factor = CGAL::square(e1.a()) + CGAL::square(e1.b());
      const FT disc = line_factor*e2.r_sq() -
                      CGAL::square(e1.a()*e2.a() + e1.b()*e2.b() + e1.c());
      CGAL::Sign sign_disc = CGAL::sign(disc);

      if (sign_disc == NEGATIVE) return res;

      const FT aux = e1.b()*e2.a() - e1.a()*e2.b();
      const FT x_base = (aux*e1.b() - e1.a()*e1.c()) / line_factor;
      const FT y_base = (-aux*e1.a() - e1.b()*e1.c()) / line_factor;

      if (sign_disc == ZERO) {
        *res++ = std::make_pair
	  ( Root_for_circles_2_2(Root_of_2(x_base), 
                                 Root_of_2(y_base)), 2u);
	return res;
      }

      // We have two intersection points, whose coordinates are one-root numbers. 
      const FT x_root_coeff = e1.b() / line_factor;
      const FT y_root_coeff = e1.a() / line_factor;

      if (CGAL::sign(e1.b()) == POSITIVE) {
        *res++ = std::make_pair
	( Root_for_circles_2_2(make_root_of_2(x_base, -x_root_coeff, disc), 
                               make_root_of_2(y_base, y_root_coeff, disc)), 1u);
        *res++ = std::make_pair
	( Root_for_circles_2_2(make_root_of_2(x_base, x_root_coeff, disc), 
                               make_root_of_2(y_base, -y_root_coeff, disc)), 1u);
      } else {
        *res++ = std::make_pair
	( Root_for_circles_2_2(make_root_of_2(x_base, x_root_coeff, disc), 
                               make_root_of_2(y_base, -y_root_coeff, disc)), 1u);
        *res++ = std::make_pair
	( Root_for_circles_2_2(make_root_of_2(x_base, -x_root_coeff, disc), 
                               make_root_of_2(y_base, y_root_coeff, disc)), 1u);
      }
      return res;
    }
  }

  template < class AK, class OutputIterator >
  inline 
  OutputIterator
  solve( const typename AK::Polynomial_for_circles_2_2 & e1,
	 const typename AK::Polynomial_1_2 & e2,
	 OutputIterator res )
  {
    return solve<AK> (e2, e1, res);
  }

  template < class AK, class OutputIterator >
  inline 
  OutputIterator
  solve( const typename AK::Polynomial_1_2 & e1,
	 const typename AK::Polynomial_1_2 & e2,
	 OutputIterator res )
  {
    typedef typename AK::FT FT;
    typedef typename AK::Root_of_2 Root_of_2;
    typedef typename AK::Root_for_circles_2_2 Root_for_circles_2_2;
    //parallele case
    const FT delta = e1.a()*e2.b() - e2.a()*e1.b();
    if(is_zero(delta)) return res;
    //case : e2 horizontal
    if(is_zero(e2.a())){
      const FT sol = -e2.c()/e2.b();
      *res++ = std::make_pair
	  ( Root_for_circles_2_2(Root_of_2(-(e1.b()*sol + e1.c())/e1.a()),
                                 Root_of_2(sol)), 1u);
      return res;
    }
    //general case
    const FT sol = (e2.a()*e1.c() - e2.c()*e1.a()) / delta; 
    *res++ = std::make_pair
	  ( Root_for_circles_2_2(Root_of_2(-(e2.b()*sol + e2.c())/e2.a()),
                                 Root_of_2(sol)), 1u);
    return res;
  }

    template < class AK >
    inline 
    Sign sign_at( const typename AK::Polynomial_1_2 & equation,
		  const typename AK::Root_for_circles_2_2 & r)
    {
      Comparison_result c = compare(r.x()*equation.a(), 
                                    -equation.c() - r.y()*equation.b());
      if(c == EQUAL) return ZERO;
      if(c == LARGER) return POSITIVE;
      return NEGATIVE;
    }


  
  } // namespace AlgebraicFunctors
} // namespace CGAL

#endif //  CGAL_ALGEBRAIC_KERNEL_FOR_CIRCLES_FUNCTIONS_ON_ROOTS_AND_POLYNOMIAL_1_2_AND_2_2_H
