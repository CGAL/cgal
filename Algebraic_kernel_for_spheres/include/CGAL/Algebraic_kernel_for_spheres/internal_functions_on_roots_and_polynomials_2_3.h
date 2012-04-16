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
//
// Author(s) : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//             Sylvain Pion
//             Pedro Machado

#ifndef CGAL_ALGEBRAIC_KERNEL_FUNCTIONS_ON_ROOTS_AND_POLYNOMIALS_2_3_H
#define CGAL_ALGEBRAIC_KERNEL_FUNCTIONS_ON_ROOTS_AND_POLYNOMIALS_2_3_H

#include <CGAL/Algebraic_kernel_for_spheres/internal_functions_on_roots_and_polynomial_1_3_and_2_3.h>

namespace CGAL {
  namespace AlgebraicSphereFunctors {

  template < class AK, class OutputIterator >
  inline 
  OutputIterator
    solve( const typename AK::Polynomial_for_spheres_2_3 &e1,
           const typename AK::Polynomial_for_spheres_2_3 &e2,
	   const typename AK::Polynomial_for_spheres_2_3 &e3,
	   OutputIterator res )
  {
	  typedef typename AK::FT FT;
    CGAL_kernel_precondition(!((e1 == e2) && (e2 == e3)));
    // we put as a precondition that the polynomial for spheres represents
    // a sphere and not an isolated point or an empty_space
    CGAL_kernel_precondition(!(e1.empty_space() || e1.isolated_point())); 
    CGAL_kernel_precondition(!(e2.empty_space() || e2.isolated_point())); 
    CGAL_kernel_precondition(!(e3.empty_space() || e3.isolated_point())); 
    typedef typename AK::Polynomial_1_3 Polynomial_1_3;
    // The degenerated cases are 2 tangent spheres
    // os 2 non-intersecting spheres
    // beacause we cannot have infinitely many solutions
    if(e1 == e2) {
      if(tangent<AK>(e1,e3)) {
        Polynomial_1_3 p = plane_from_2_spheres<AK>(e1,e3);
        return internal::solve_tangent<AK>(p,e1,res);
      }
      CGAL_kernel_precondition(!(intersect<AK>(e1,e3)));
      return res;
    }
    if((e1 == e3) || (e2 == e3)) {
      if(tangent<AK>(e1,e2)) {
        Polynomial_1_3 p = plane_from_2_spheres<AK>(e1,e2);
        return internal::solve_tangent<AK>(p,e1,res);
      }
      CGAL_kernel_precondition(!(intersect<AK>(e1,e2)));
      return res;
    }
    
    // non degenerated case
    if(intersect<AK>(e1,e2)) {
      Polynomial_1_3 p1 = plane_from_2_spheres<AK>(e1,e2);
      if(intersect<AK>(e2,e3)) {
        Polynomial_1_3 p2 = plane_from_2_spheres<AK>(e2,e3); 
        if(same_solutions<FT>(p1,p2)) {
	        const FT sq_d1 = CGAL::square(p1.a()*e1.a() + p1.b()*e1.b() +
	                                p1.c()*e1.c() + p1.d()) /
							(square(p1.a()) + square(p1.b()) + square(p1.c()));
					const FT r1_sqr = e1.r_sq() - sq_d1;

	        const FT sq_d2 = CGAL::square(p2.a()*e2.a() + p2.b()*e2.b() +
	                                p2.c()*e2.c() + p2.d()) /
							(square(p2.a()) + square(p2.b()) + square(p2.c()));
					const FT r2_sqr = e2.r_sq() - sq_d2;
					if(r1_sqr != r2_sqr) return res;
	        // otherwise there are an infinite number of points
	        // this is not allowed
					CGAL_kernel_precondition(r1_sqr == 0);
					return internal::solve_tangent<AK>(p1,e1,res);
	      }        
				return solve<AK>(p1,p2,e2,res);
      } return res;
    } return res;
  }

  template <class AK>
  typename AK::Root_for_spheres_2_3
  x_critical_point(const typename AK::Polynomial_for_spheres_2_3 & s, bool i)
  {
    typedef typename AK::Root_of_2 Root_of_2;
    typedef typename AK::Root_for_spheres_2_3 Root_for_spheres_2_3; 

    return Root_for_spheres_2_3(
            make_root_of_2(s.a(),typename AK::FT(i?-1:1),s.r_sq()),
            Root_of_2(s.b()),
            Root_of_2(s.c()));
  }

  template <class AK, class OutputIterator>
  OutputIterator
  x_critical_points(const typename AK::Polynomial_for_spheres_2_3 & s, OutputIterator res)
  {
    typedef typename AK::Root_of_2 Root_of_2;
    typedef typename AK::Root_for_spheres_2_3 Root_for_spheres_2_3; 
    typedef typename AK::FT FT;
    
    *res++ =  Root_for_spheres_2_3(make_root_of_2(s.a(),FT(-1),s.r_sq()),
                                Root_of_2(s.b()),
                                Root_of_2(s.c()));
    *res++ =  Root_for_spheres_2_3(make_root_of_2(s.a(),FT(1),s.r_sq()),
                                Root_of_2(s.b()),
                                Root_of_2(s.c()));
    return res;
  }

  template <class AK>
  typename AK::Root_for_spheres_2_3
  y_critical_point(const typename AK::Polynomial_for_spheres_2_3 &s, bool i)
  {
    typedef typename AK::Root_of_2 Root_of_2;
    typedef typename AK::Root_for_spheres_2_3 Root_for_spheres_2_3; 

    return Root_for_spheres_2_3(
            Root_of_2(s.a()),
            make_root_of_2(s.b(),typename AK::FT(i?-1:1),s.r_sq()),
            Root_of_2(s.c()));
  }
  
  template <class AK, class OutputIterator>
  OutputIterator
  y_critical_points(const typename AK::Polynomial_for_spheres_2_3 & s, OutputIterator res)
  {
    typedef typename AK::Root_of_2 Root_of_2;
    typedef typename AK::Root_for_spheres_2_3 Root_for_spheres_2_3;
    typedef typename AK::FT FT;

    *res++ =  Root_for_spheres_2_3(Root_of_2(s.a()),
                                make_root_of_2(s.b(),FT(-1),s.r_sq()),
                                Root_of_2(s.c()));
    *res++ =  Root_for_spheres_2_3(Root_of_2(s.a()),
                                make_root_of_2(s.b(),FT(1),s.r_sq()),
                                Root_of_2(s.c()));
    return res;
  }

  template <class AK>
  typename AK::Root_for_spheres_2_3
  z_critical_point(const typename AK::Polynomial_for_spheres_2_3 &s, bool i)
  {
    typedef typename AK::Root_of_2 Root_of_2;
    typedef typename AK::Root_for_spheres_2_3 Root_for_spheres_2_3; 

    return Root_for_spheres_2_3(Root_of_2(s.a()),
                                Root_of_2(s.b()),
            make_root_of_2(s.c(),typename AK::FT(i?-1:1),s.r_sq()));
  }
  
  template <class AK, class OutputIterator>
  OutputIterator
  z_critical_points(const typename AK::Polynomial_for_spheres_2_3 & s, OutputIterator res)
  {
    typedef typename AK::Root_of_2 Root_of_2;
    typedef typename AK::Root_for_spheres_2_3 Root_for_spheres_2_3;
    typedef typename AK::FT FT;

    *res++ =  Root_for_spheres_2_3(Root_of_2(s.a()),
                                Root_of_2(s.b()),
                                make_root_of_2(s.c(),FT(-1),s.r_sq()));
    *res++ =  Root_for_spheres_2_3(Root_of_2(s.a()),
                                Root_of_2(s.b()),
                                make_root_of_2(s.c(),FT(1),s.r_sq()));
    return res;
  }

  template <class AK>
  typename AK::Root_for_spheres_2_3
  x_critical_point( const std::pair<typename AK::Polynomial_for_spheres_2_3, 
                                     typename AK::Polynomial_1_3 > &c, bool i)
  {
    typedef typename AK::FT FT;
    typedef typename AK::Root_of_2 Root_of_2;
    typedef typename AK::Root_for_spheres_2_3 Root_for_spheres_2_3; 
    typedef typename AK::Polynomial_for_spheres_2_3 Polynomial_for_spheres_2_3;
    typedef typename AK::Polynomial_1_3 Polynomial_1_3;

    const Polynomial_for_spheres_2_3 &s = c.first;
    const Polynomial_1_3 &p = c.second;

    // It has to be the equation of a diametral circle
    CGAL_kernel_precondition((intersect<AK>(p,s)));
    CGAL_kernel_precondition(CGAL_NTS sign(p.a() * s.a() + p.b() * s.b() + 
                                      p.c() * s.c() + p.d()) == ZERO);
    CGAL_kernel_precondition(!(is_zero(p.b()) && is_zero(p.c())));

    const FT sqbc = CGAL::square(p.b()) + CGAL::square(p.c());
    const FT sq_sum = sqbc + CGAL::square(p.a());
    const FT delta = (sqbc * s.r_sq())/sq_sum;

    const FT cy = (p.a()*p.b())/sqbc;
    const FT cz = (p.a()*p.c())/sqbc;

    const Root_of_2 x = make_root_of_2(s.a(),FT(i?-1:1),delta);
    const Root_of_2 y = make_root_of_2(s.b(),FT(i?(cy):(-cy)),delta);
    const Root_of_2 z = make_root_of_2(s.c(),FT(i?(cz):(-cz)),delta);

    return Root_for_spheres_2_3(x,y,z);
  }

  template <class AK, class OutputIterator>
  OutputIterator
  x_critical_points( const std::pair<typename AK::Polynomial_for_spheres_2_3, 
                                     typename AK::Polynomial_1_3 > &c, 
                     OutputIterator res)
  {
    typedef typename AK::FT FT;
    typedef typename AK::Root_of_2 Root_of_2;
    typedef typename AK::Root_for_spheres_2_3 Root_for_spheres_2_3; 
    typedef typename AK::Polynomial_for_spheres_2_3 Polynomial_for_spheres_2_3;
    typedef typename AK::Polynomial_1_3 Polynomial_1_3;

    const Polynomial_for_spheres_2_3 &s = c.first;
    const Polynomial_1_3 &p = c.second;

    // It has to be the equation of a diametral circle
    CGAL_kernel_precondition((intersect<AK>(p,s)));
    CGAL_kernel_precondition(CGAL_NTS sign(p.a() * s.a() + p.b() * s.b() + 
                                           p.c() * s.c() + p.d()) == ZERO);
    CGAL_kernel_precondition(!(is_zero(p.b()) && is_zero(p.c())));

    const FT sqbc = CGAL::square(p.b()) + CGAL::square(p.c());
    const FT sq_sum = sqbc + CGAL::square(p.a());
    const FT delta = (sqbc * s.r_sq())/sq_sum;

    const FT cy = (p.a()*p.b())/sqbc;
    const FT cz = (p.a()*p.c())/sqbc;

    const Root_of_2 x1 = make_root_of_2(s.a(),-1,delta);
    const Root_of_2 y1 = make_root_of_2(s.b(),cy,delta);
    const Root_of_2 z1 = make_root_of_2(s.c(),cz,delta);
    const Root_of_2 x2 = make_root_of_2(s.a(),1,delta);
    const Root_of_2 y2 = make_root_of_2(s.b(),-cy,delta);
    const Root_of_2 z2 = make_root_of_2(s.c(),-cz,delta);

    *res++ =  Root_for_spheres_2_3(x1,y1,z1);
    *res++ =  Root_for_spheres_2_3(x2,y2,z2);
    return res;
  }

  template <class AK>
  typename AK::Root_for_spheres_2_3
  y_critical_point( const std::pair<typename AK::Polynomial_for_spheres_2_3, 
                                     typename AK::Polynomial_1_3 > &c, bool i)
  {
    typedef typename AK::FT FT;
    typedef typename AK::Root_of_2 Root_of_2;
    typedef typename AK::Root_for_spheres_2_3 Root_for_spheres_2_3; 
    typedef typename AK::Polynomial_for_spheres_2_3 Polynomial_for_spheres_2_3;
    typedef typename AK::Polynomial_1_3 Polynomial_1_3;

    const Polynomial_for_spheres_2_3 &s = c.first;
    const Polynomial_1_3 &p = c.second;

    // It has to be the equation of a diametral circle
    CGAL_kernel_precondition((intersect<AK>(p,s)));
    CGAL_kernel_precondition(CGAL_NTS sign(p.a() * s.a() + p.b() * s.b() + 
                                           p.c() * s.c() + p.d()) == ZERO);
    CGAL_kernel_precondition(!(is_zero(p.a()) && is_zero(p.c())));

    const FT sqac = CGAL::square(p.a()) + CGAL::square(p.c());
    const FT sq_sum = sqac + CGAL::square(p.b());
    const FT delta = (sqac * s.r_sq())/sq_sum;

    const FT cx = (p.a()*p.b())/sqac;
    const FT cz = (p.c()*p.b())/sqac;

    if(!is_positive(cx)) {
        const Root_of_2 x = make_root_of_2(s.a(),FT(i?(cx):(-cx)),delta);
        const Root_of_2 y = make_root_of_2(s.b(),FT(i?-1:1),delta);
        const Root_of_2 z = make_root_of_2(s.c(),FT(i?(cz):(-cz)),delta);
      return Root_for_spheres_2_3(x,y,z);
    } else {
        const Root_of_2 x = make_root_of_2(s.a(),FT(i?(-cx):(cx)),delta);
        const Root_of_2 y = make_root_of_2(s.b(),FT(i?1:-1),delta);
        const Root_of_2 z = make_root_of_2(s.c(),FT(i?(-cz):(cz)),delta);
      return Root_for_spheres_2_3(x,y,z);
    } 
  }

  template <class AK, class OutputIterator>
  OutputIterator
  y_critical_points( const std::pair<typename AK::Polynomial_for_spheres_2_3, 
                                     typename AK::Polynomial_1_3 > &c, 
                     OutputIterator res)
  {
    typedef typename AK::FT FT;
    typedef typename AK::Root_of_2 Root_of_2;
    typedef typename AK::Root_for_spheres_2_3 Root_for_spheres_2_3; 
    typedef typename AK::Polynomial_for_spheres_2_3 Polynomial_for_spheres_2_3;
    typedef typename AK::Polynomial_1_3 Polynomial_1_3;

    const Polynomial_for_spheres_2_3 &s = c.first;
    const Polynomial_1_3 &p = c.second;

    // It has to be the equation of a diametral circle
    CGAL_kernel_precondition((intersect<AK>(p,s)));
    CGAL_kernel_precondition(CGAL_NTS sign(p.a() * s.a() + p.b() * s.b() + 
                                           p.c() * s.c() + p.d()) == ZERO);
    CGAL_kernel_precondition(!(is_zero(p.a()) && is_zero(p.c())));

    const FT sqac = CGAL::square(p.a()) + CGAL::square(p.c());
    const FT sq_sum = sqac + CGAL::square(p.b());
    const FT delta = (sqac * s.r_sq())/sq_sum;

    const FT cx = (p.a()*p.b())/sqac;
    const FT cz = (p.c()*p.b())/sqac;

    const Root_of_2 x1 = make_root_of_2(s.a(),cx,delta);
    const Root_of_2 y1 = make_root_of_2(s.b(),FT(-1),delta);
    const Root_of_2 z1 = make_root_of_2(s.c(),cz,delta);
    const Root_of_2 x2 = make_root_of_2(s.a(),-cx,delta);
    const Root_of_2 y2 = make_root_of_2(s.b(),FT(1),delta);
    const Root_of_2 z2 = make_root_of_2(s.c(),-cz,delta);

    if(!is_positive(cx)) {
      *res++ =  Root_for_spheres_2_3(x1,y1,z1);
      *res++ =  Root_for_spheres_2_3(x2,y2,z2);
    } else {
      *res++ =  Root_for_spheres_2_3(x2,y2,z2);
      *res++ =  Root_for_spheres_2_3(x1,y1,z1);
    } 
    return res;
  }

  template <class AK>
  typename AK::Root_for_spheres_2_3
  z_critical_point( const std::pair<typename AK::Polynomial_for_spheres_2_3, 
                                     typename AK::Polynomial_1_3 > &c, bool i)
  {
    typedef typename AK::FT FT;
    typedef typename AK::Root_of_2 Root_of_2;
    typedef typename AK::Root_for_spheres_2_3 Root_for_spheres_2_3; 
    typedef typename AK::Polynomial_for_spheres_2_3 Polynomial_for_spheres_2_3;
    typedef typename AK::Polynomial_1_3 Polynomial_1_3;

    const Polynomial_for_spheres_2_3 &s = c.first;
    const Polynomial_1_3 &p = c.second;

    // It has to be the equation of a diametral circle
    CGAL_kernel_precondition((intersect<AK>(p,s)));
    CGAL_kernel_precondition(CGAL_NTS sign(p.a() * s.a() + p.b() * s.b() + 
                                           p.c() * s.c() + p.d()) == ZERO);
    CGAL_kernel_precondition(!(is_zero(p.a()) && is_zero(p.b())));

    const FT sqab = CGAL::square(p.a()) + CGAL::square(p.b());
    const FT sq_sum = sqab + CGAL::square(p.c());
    const FT delta = (sqab * s.r_sq())/sq_sum;

    const FT cx = (p.a()*p.c())/sqab;
    const FT cy = (p.c()*p.b())/sqab;

    if(is_negative(cx)) {
        const Root_of_2 x = make_root_of_2(s.a(),FT(i?(cx):(-cx)),delta);
        const Root_of_2 y = make_root_of_2(s.b(),FT(i?(cy):(-cy)),delta);
        const Root_of_2 z = make_root_of_2(s.c(),FT(i?-1:1),delta);
      return Root_for_spheres_2_3(x,y,z);
    } else if(is_zero(cx)) {
      if(!is_positive(cy)) {
        const Root_of_2 x = s.a();
        const Root_of_2 y = make_root_of_2(s.b(),FT(i?(cy):(-cy)),delta);
        const Root_of_2 z = make_root_of_2(s.c(),FT(i?-1:1),delta);
        return Root_for_spheres_2_3(x,y,z);
      } else {
        const Root_of_2 x = s.a();
        const Root_of_2 y = make_root_of_2(s.b(),FT(i?(-cy):(cy)),delta);
        const Root_of_2 z = make_root_of_2(s.c(),FT(i?1:-1),delta);
        return Root_for_spheres_2_3(x,y,z);
      }
    } else {
        const Root_of_2 x = make_root_of_2(s.a(),FT(i?(-cx):(cx)),delta);
        const Root_of_2 y = make_root_of_2(s.b(),FT(i?(-cy):(cy)),delta);
        const Root_of_2 z = make_root_of_2(s.c(),FT(i?1:-1),delta);
      return Root_for_spheres_2_3(x,y,z);
    } 
  }

  template <class AK, class OutputIterator>
  OutputIterator
  z_critical_points( const std::pair<typename AK::Polynomial_for_spheres_2_3, 
                                     typename AK::Polynomial_1_3 > &c, 
                     OutputIterator res)
  {
    typedef typename AK::FT FT;
    typedef typename AK::Root_of_2 Root_of_2;
    typedef typename AK::Root_for_spheres_2_3 Root_for_spheres_2_3; 
    typedef typename AK::Polynomial_for_spheres_2_3 Polynomial_for_spheres_2_3;
    typedef typename AK::Polynomial_1_3 Polynomial_1_3;

    const Polynomial_for_spheres_2_3 &s = c.first;
    const Polynomial_1_3 &p = c.second;

    // It has to be the equation of a diametral circle
    CGAL_kernel_precondition((intersect<AK>(p,s)));
    CGAL_kernel_precondition(CGAL_NTS sign(p.a() * s.a() + p.b() * s.b() + 
                                           p.c() * s.c() + p.d()) == ZERO);
    CGAL_kernel_precondition(!(is_zero(p.a()) && is_zero(p.b())));

    const FT sqab = CGAL::square(p.a()) + CGAL::square(p.b());
    const FT sq_sum = sqab + CGAL::square(p.c());
    const FT delta = (sqab * s.r_sq())/sq_sum;

    const FT cx = (p.a()*p.c())/sqab;
    const FT cy = (p.c()*p.b())/sqab;

    if(is_negative(cx)) {
      const Root_of_2 x1 = make_root_of_2(s.a(),(cx),delta);
      const Root_of_2 y1 = make_root_of_2(s.b(),(cy),delta);
      const Root_of_2 z1 = make_root_of_2(s.c(),-1,delta);
      const Root_of_2 x2 = make_root_of_2(s.a(),(-cx),delta);
      const Root_of_2 y2 = make_root_of_2(s.b(),(-cy),delta);
      const Root_of_2 z2 = make_root_of_2(s.c(),1,delta);
      *res++ = Root_for_spheres_2_3(x1,y1,z1);
      *res++ = Root_for_spheres_2_3(x2,y2,z2);
    } else if(is_zero(cx)) {
      if(!is_positive(cy)) {
        const Root_of_2 x1 = s.a();
        const Root_of_2 y1 = make_root_of_2(s.b(),(cy),delta);
        const Root_of_2 z1 = make_root_of_2(s.c(),FT(-1),delta);
        const Root_of_2 y2 = make_root_of_2(s.b(),(-cy),delta);
        const Root_of_2 z2 = make_root_of_2(s.c(),FT(1),delta);
        *res++ = Root_for_spheres_2_3(x1,y1,z1);
        *res++ = Root_for_spheres_2_3(x1,y2,z2);
      } else {
        const Root_of_2 x1 = s.a();
        const Root_of_2 y1 = make_root_of_2(s.b(),(-cy),delta);
        const Root_of_2 z1 = make_root_of_2(s.c(),FT(1),delta);
        const Root_of_2 y2 = make_root_of_2(s.b(),(cy),delta);
        const Root_of_2 z2 = make_root_of_2(s.c(),FT(-1),delta);
        *res++ = Root_for_spheres_2_3(x1,y1,z1);
        *res++ = Root_for_spheres_2_3(x1,y2,z2);
      }
    } else {
      const Root_of_2 x1 = make_root_of_2(s.a(),(-cx),delta);
      const Root_of_2 y1 = make_root_of_2(s.b(),(-cy),delta);
      const Root_of_2 z1 = make_root_of_2(s.c(),FT(1),delta);
      const Root_of_2 x2 = make_root_of_2(s.a(),(cx),delta);
      const Root_of_2 y2 = make_root_of_2(s.b(),(cy),delta);
      const Root_of_2 z2 = make_root_of_2(s.c(),FT(-1),delta);
      *res++ = Root_for_spheres_2_3(x1,y1,z1);
      *res++ = Root_for_spheres_2_3(x2,y2,z2);
    }
    return res;
  }

  } // namespace AlgebraicSphereFunctors
} // namespace CGAL

#endif //  CGAL_ALGEBRAIC_KERNEL_FUNCTIONS_ON_ROOTS_AND_POLYNOMIALS_2_3_H
