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
// Author(s)     : Monique Teillaud, Sylvain Pion, Julien Hazebrouck

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_ALGEBRAIC_KERNEL_FUNCTIONS_ON_ROOTS_AND_POLYNOMIAL_1_2_AND_2_2_H
#define CGAL_ALGEBRAIC_KERNEL_FUNCTIONS_ON_ROOTS_AND_POLYNOMIAL_1_2_AND_2_2_H


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
      
      const FT a = 1;
      const FT b = -2*e2.a();
      const FT c_t = e1.c()/e1.b();
      const FT c = CGAL::square(e2.a()) +
	CGAL::square(e2.b() + c_t)
	- e2.r_sq();
      const FT cond = CGAL::square(b) - 4 * c;
    
      if(cond < 0)
	return res;
      if (cond == 0) {	
	*res++ = std::make_pair
	  ( Root_for_circles_2_2(Root_of_2(e2.a()), 
                                 Root_of_2(-c_t)), 2u);
	return res;
      }
      const Root_of_2 x_res1 = make_root_of_2(a, b, c, true, true);
      const Root_of_2 x_res2 = make_root_of_2(a, b, c, false, true);
      const Root_of_2 y_res = Root_of_2(-c_t); 
      *res++ = std::make_pair
	( Root_for_circles_2_2(x_res1, y_res), 1u);
      *res++ = std::make_pair
	( Root_for_circles_2_2(x_res2,y_res), 1u);
      return res;
    }
    else if(is_zero(e1.b())){//vertical line
       
      const FT a = 1;
      const FT b = -2*e2.b();
      const FT c_t = e1.c()/e1.a();
      const FT c = CGAL::square(e2.b()) +
	CGAL::square(c_t + e2.a())
	- e2.r_sq();
      const FT cond = CGAL::square(b) - 4 * c;
    
      if(cond < 0)
	return res;
      if (cond == 0) {
	*res++ = std::make_pair
	  ( Root_for_circles_2_2(Root_of_2(-c_t), 
                                 Root_of_2(e2.b())), 2u);
	return res;
      }
      const Root_of_2 y_res1 = make_root_of_2(a, b, c, true, true);
      const Root_of_2 y_res2 = make_root_of_2(a, b, c, false, true);
      Root_of_2 x_res = Root_of_2(-c_t); 
      *res++ = std::make_pair
	( Root_for_circles_2_2(x_res, y_res1), 1u);
      *res++ = std::make_pair
	( Root_for_circles_2_2(x_res,y_res2), 1u);
      return res;
    }
    else{
      //general case
      const FT a1_sq = CGAL::square(e1.a());
      const FT b1_sq = CGAL::square(e1.b());
      const FT b2_sq = CGAL::square(e2.b());
      
      const FT a = b1_sq + a1_sq;
      const FT b = 2*e1.c()*e1.b() + 2*e2.a()*e1.a()*e1.b() - 2*a1_sq*e2.b();
      const FT c = CGAL::square(e2.a()*e1.a() + e1.c()) + a1_sq*b2_sq - a1_sq*e2.r_sq();
      
      const FT cond = CGAL::square(b) - 4 *a*c;
    
      if(cond < 0)
	return res;


      if (cond == 0) {
        const FT sol = -b/(2*a);
	*res++ = std::make_pair
	  ( Root_for_circles_2_2(Root_of_2((e1.b()*sol + e1.c()) / -e1.a()), 
                                 Root_of_2(sol)), 2u);
	return res;
      }
    
      const Root_of_2 y_res1 = make_root_of_2(a, b, c, true, true);
      const Root_of_2 x_res1 = (-e1.b()/e1.a()) * y_res1 - (e1.c()/e1.a());
      const Root_of_2 y_res2 = make_root_of_2(a, b, c, false, true);
      const Root_of_2 x_res2 = (-e1.b()/e1.a()) * y_res2 - (e1.c()/e1.a());

      if(x_res2 < x_res1){
	*res++ = std::make_pair
	  ( Root_for_circles_2_2(x_res2, y_res2), 1u);
	*res++ = std::make_pair
	  ( Root_for_circles_2_2(x_res1,y_res1), 1u);
	return res;
      }
      else{
        *res++ = std::make_pair
	  ( Root_for_circles_2_2(x_res1, y_res1), 1u);
	*res++ = std::make_pair
	  ( Root_for_circles_2_2(x_res2,y_res2), 1u);
	return res;
      }     
    }
  }

  template < class AK, class OutputIterator >
  inline 
  OutputIterator
  solve( const typename AK::Polynomial_for_circles_2_2 & e1,
	 const typename AK::Polynomial_1_2 & e2,
	 OutputIterator res )
  {
    return solve<AK> ( e2, e1, res);
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
      typedef typename AK::Root_of_2 Root_of_2;
      Comparison_result c = compare(r.x()*equation.a(), 
                                    -equation.c() - r.y()*equation.b());
      if(c == EQUAL) return ZERO;
      if(c == LARGER) return POSITIVE;
      return NEGATIVE;
    }


  
  } // namespace AlgebraicFunctors
} // namespace CGAL

#endif //  CGAL_ALGEBRAIC_KERNEL_FUNCTIONS_ON_ROOTS_AND_POLYNOMIAL1_2_AND_2_2_H
