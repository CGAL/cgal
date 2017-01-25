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

#ifndef CGAL_ALGEBRAIC_KERNEL_FUNCTIONS_ON_ROOTS_AND_POLYNOMIAL_1_3_AND_2_3_H
#define CGAL_ALGEBRAIC_KERNEL_FUNCTIONS_ON_ROOTS_AND_POLYNOMIAL_1_3_AND_2_3_H

#include <CGAL/license/Circular_kernel_3.h>


#include <vector>
#include <iterator> // for std::back_inserter
#include <utility>  // for std::pair and std::make_pair
#include <CGAL/enum.h> // for CGAL::Sign
#include <CGAL/use.h>  // for CGAL_USE_TYPE()

namespace CGAL {
  namespace AlgebraicSphereFunctors {

  template < class AK >
  inline 
  Sign sign_at( const typename AK::Polynomial_for_spheres_2_3 & equation,
		  const typename AK::Root_for_spheres_2_3 &r){
    return CGAL_NTS sign(r.evaluate(equation));
  }

  template < class AK >
  inline 
  Sign sign_at(const typename AK::Polynomial_1_3 & equation,
		const typename AK::Root_for_spheres_2_3 &r){
    return CGAL_NTS sign(r.evaluate(equation));
  }

  template < class AK >
  inline 
  bool contains(const typename AK::Polynomials_for_line_3 & equation,
		const typename AK::Root_for_spheres_2_3 &r){
    return r.is_on_line(equation); 
  }

  template <class AK>
  bool intersect(const typename AK::Polynomial_1_3 & p1, 
                 const typename AK::Polynomial_1_3 & p2) {
    typedef typename AK::RT RT;

    CGAL_USE_TYPE(RT);
    CGAL_kernel_precondition(!(same_solutions<RT>(p1,p2)));

    if(p1.empty_space()) return false;
    if(p2.empty_space()) return false;

    return !((p2.b() * p1.a() == p1.b() * p2.a()) &&
             (p2.c() * p1.b() == p1.c() * p2.b()) &&
             (p2.c() * p1.a() == p1.c() * p2.a()));
  }

  template <class AK>
  inline
  typename AK::Polynomials_for_line_3
  line_from_2_planes(const typename AK::Polynomial_1_3 & p1, 
                     const typename AK::Polynomial_1_3 & p2) 
  {
    typedef typename AK::Polynomials_for_line_3 Polynomials_for_line_3;
    typedef typename AK::FT FT;
    CGAL_kernel_precondition(intersect<AK>(p1,p2));

    const FT a1 = p1.b() * p2.c() - p1.c() * p2.b();
    const FT a2 = p1.c() * p2.a() - p1.a() * p2.c();
    const FT a3 = p1.a() * p2.b() - p1.b() * p2.a();

    if(!is_zero(a1)) {
      const FT b1 = 0;
      const FT b2 = (p2.d() * p1.c() - p2.c() * p1.d()) / a1;
      const FT b3 = (p2.b() * p1.d() - p2.d() * p1.b()) / a1;
      return Polynomials_for_line_3(a1, b1, a2, b2, a3, b3);
    } 

    if(!is_zero(a2)) {
      const FT b1 = (p1.d() * p2.c() - p1.c() * p2.d()) / a2;
      const FT b2 = 0;
      const FT b3 = (p1.a() * p2.d() - p1.d() * p2.a()) / a2;
      return Polynomials_for_line_3(a1, b1, a2, b2, a3, b3);
    }
    // a3 must not be 0
    CGAL_kernel_precondition(!is_zero(a3));

    const FT b1 = (p2.d() * p1.b() - p2.b() * p1.d()) / a3;
    const FT b2 = (p2.a() * p1.d() - p2.d() * p1.a()) / a3;
    const FT b3 = 0;
    return Polynomials_for_line_3(a1, b1, a2, b2, a3, b3);
 
  }

  template <class AK>
  inline
  bool intersect(const typename AK::Polynomial_for_spheres_2_3 & s1, 
                 const typename AK::Polynomial_for_spheres_2_3 & s2) {
    typedef typename AK::FT FT;
    typedef typename AK::Root_of_2 Root_of_2;
    
    const FT dx = s2.a() - s1.a();
    const FT dy = s2.b() - s1.b();
    const FT dz = s2.c() - s1.c();
    const FT d2 = CGAL::square(dx) +
                  CGAL::square(dy) +
                  CGAL::square(dz);
    const FT sq_r1r2 = s1.r_sq()*s2.r_sq();
    const FT sq_r1_p_sq_r2 = s1.r_sq() + s2.r_sq();
  
    Root_of_2 left_1 = make_root_of_2(d2,FT(-2),sq_r1r2);
    if(left_1 > sq_r1_p_sq_r2) return false;
    Root_of_2 left_2 = make_root_of_2(d2,FT(2),sq_r1r2);
    if(left_2 < sq_r1_p_sq_r2) return false;
    return true;
  }

  template <class AK>
  inline
  bool intersect(const typename AK::Polynomial_1_3 & p, 
                 const typename AK::Polynomial_for_spheres_2_3 & s) {
    return CGAL::square(p.a()*s.a() + p.b()*s.b() + p.c()*s.c() + p.d()) <=  
        ((CGAL::square(p.a()) + CGAL::square(p.b()) + CGAL::square(p.c())) * 
        s.r_sq());
  }

  template <class AK>
  inline
  bool tangent(const typename AK::Polynomial_for_spheres_2_3 & s1, 
                 const typename AK::Polynomial_for_spheres_2_3 & s2) {
    typedef typename AK::RT RT;
    typedef typename AK::Root_of_2 Root_of_2;
    
    const RT dx = s2.a() - s1.a();
    const RT dy = s2.b() - s1.b();
    const RT dz = s2.c() - s1.c();
    const RT d2 = CGAL::square(dx) +
                  CGAL::square(dy) +
                  CGAL::square(dz);
    const RT sq_r1r2 = s1.r_sq()*s2.r_sq();
    const RT sq_r1_p_sq_r2 = s1.r_sq() + s2.r_sq();
  
    Root_of_2 left_1 = make_root_of_2(d2,RT(-2),sq_r1r2);
    if(left_1 == sq_r1_p_sq_r2) return true;
    Root_of_2 left_2 = make_root_of_2(d2,RT(2),sq_r1r2);
    if(left_2 == sq_r1_p_sq_r2) return true;
    return false;
  }

  template <class AK>
  inline
  bool tangent(const typename AK::Polynomial_1_3  & p, 
                 const typename AK::Polynomial_for_spheres_2_3 & s) {
    return CGAL::square(p.a()*s.a() + p.b()*s.b() + p.c()*s.c() + p.d()) ==  
        ((CGAL::square(p.a()) + CGAL::square(p.b()) + CGAL::square(p.c())) * 
        s.r_sq());
  }

  template <class AK>
  typename AK::Polynomial_1_3
  plane_from_2_spheres(const typename AK::Polynomial_for_spheres_2_3 & s1, 
                       const typename AK::Polynomial_for_spheres_2_3 & s2) 
  {
    typedef typename AK::Polynomial_1_3 Polynomial_1_3;
    typedef typename AK::RT RT;

    // In order to have a radical plane
    // the spheres given must intersect
    CGAL_kernel_precondition(intersect<AK>(s1,s2));

    const RT a = 2*(s2.a() - s1.a());
    const RT b = 2*(s2.b() - s1.b());
    const RT c = 2*(s2.c() - s1.c());
    const RT d = CGAL::square(s1.a()) + 
           CGAL::square(s1.b()) +
           CGAL::square(s1.c()) - s1.r_sq() - 
      CGAL::square(s2.a()) - 
      CGAL::square(s2.b()) -
      CGAL::square(s2.c()) + s2.r_sq();
    return Polynomial_1_3(a, b, c, d);
  }

  template < class AK, class OutputIterator >
  inline 
  OutputIterator
    solve(const typename AK::Polynomials_for_line_3 &p,
          const typename AK::Polynomial_for_spheres_2_3& s,
	 OutputIterator res )
  {
    typedef typename AK::FT FT;
    typedef typename AK::Root_for_spheres_2_3 Root_for_spheres_2_3;
    typedef typename AK::Root_of_2 Root_of_2;

    // we must have a line
    CGAL_kernel_precondition(!p.degenerated());  

    const FT sq_a1 = CGAL::square(p.a1());
    const FT sq_a2 = CGAL::square(p.a2());
    const FT sq_a3 = CGAL::square(p.a3());
    const FT sq_a = CGAL::square(s.a());
    const FT sq_b = CGAL::square(s.b());
    const FT sq_c = CGAL::square(s.c());

    const FT a = sq_a1 + sq_a2 + sq_a3;
    const FT b = p.a1() * (p.b1() - s.a()) +
                 p.a2() * (p.b2() - s.b()) +
                 p.a3() * (p.b3() - s.c());
    const FT c = CGAL::square(p.b1()) + 
                 CGAL::square(p.b2()) + 
                 CGAL::square(p.b3()) + 
                 sq_a + sq_b + sq_c - 
                 2*(s.a() * p.b1() + 
                    s.b() * p.b2() +
                    s.c() * p.b3()) - s.r_sq();
    const FT alpha = -b/a;
    const FT gama = CGAL::square(alpha) - (c/a);

    if(is_negative(gama)) return res;

    if(is_zero(gama)) {
      *res++ = std::make_pair(
        Root_for_spheres_2_3(Root_of_2(p.a1() * alpha + p.b1()),
                             Root_of_2(p.a2() * alpha + p.b2()),
                             Root_of_2(p.a3() * alpha + p.b3())), 
        2); 
      return res;
    }

    const Root_of_2 t1 = make_root_of_2(alpha,FT(-1),gama);
    const Root_of_2 t2 = make_root_of_2(alpha,FT(1),gama);

    bool first_t1 = true;
    Sign sign_a1 = CGAL_NTS sign(p.a1());
    Sign sign_a2 = CGAL_NTS sign(p.a2());
    Sign sign_a3 = CGAL_NTS sign(p.a3());
    if(sign_a1 == ZERO) {
      if(sign_a2 == ZERO) {
        first_t1 = (sign_a3 == POSITIVE);
      } else first_t1 = (sign_a2 == POSITIVE);
    } else first_t1 = (sign_a1 == POSITIVE);

    if(first_t1) {
      *res++ = std::make_pair(
        Root_for_spheres_2_3(p.a1() * t1 + p.b1(),
                             p.a2() * t1 + p.b2(),
                             p.a3() * t1 + p.b3()), 
        1); 
      *res++ = std::make_pair(
        Root_for_spheres_2_3(p.a1() * t2 + p.b1(),
                             p.a2() * t2 + p.b2(),
                             p.a3() * t2 + p.b3()), 
        1);
    } else {
      *res++ = std::make_pair(
        Root_for_spheres_2_3(p.a1() * t2 + p.b1(),
                             p.a2() * t2 + p.b2(),
                             p.a3() * t2 + p.b3()), 
        1); 
      *res++ = std::make_pair(
        Root_for_spheres_2_3(p.a1() * t1 + p.b1(),
                             p.a2() * t1 + p.b2(),
                             p.a3() * t1 + p.b3()), 
        1);
    } 

    return res;
  }

  template < class AK, class OutputIterator >
  inline 
  OutputIterator
    solve(const typename AK::Polynomial_for_spheres_2_3& s,
          const typename AK::Polynomials_for_line_3 &p,
	 OutputIterator res )
  {
    return solve<AK>(p,s,res);
  }

  namespace internal {

    template < class AK, class OutputIterator >
    inline 
    OutputIterator
      solve_tangent(const typename AK::Polynomial_1_3 & p, 
            const typename AK::Polynomial_for_spheres_2_3& s,
  	    OutputIterator res )
    {
      typedef typename AK::FT FT;
      typedef typename AK::Root_of_2 Root_of_2;
      typedef typename AK::Root_for_spheres_2_3 Root_for_spheres_2_3;
      // It returns no solution when there is infinitely solution
      // For now, this solve is only for internal computations purpose
      CGAL_kernel_precondition(
        CGAL::square(p.a()*s.a() + p.b()*s.b() + p.c()*s.c() + p.d()) ==  
        ((CGAL::square(p.a()) + CGAL::square(p.b()) + CGAL::square(p.c())) * 
        s.r_sq()));

      const FT t = -(p.a()*s.a() + p.b()*s.b() + p.c()*s.c() + p.d()) / 
        (CGAL::square(p.a()) + CGAL::square(p.b()) + CGAL::square(p.c())); 

      *res++ = std::make_pair(
        Root_for_spheres_2_3(Root_of_2(p.a() * t + s.a()),
                             Root_of_2(p.b() * t + s.b()),
                             Root_of_2(p.c() * t + s.c())), 
        2); 

      return res;
    }
  }

  template < class AK, class OutputIterator >
  inline 
  OutputIterator
    solve(const typename AK::Polynomial_1_3 & e1,
	  const typename AK::Polynomial_1_3 & e2,
          const typename AK::Polynomial_for_spheres_2_3& e3,
	 OutputIterator res )
  {
    typedef typename AK::RT RT;
    typedef typename AK::Polynomials_for_line_3 Polynomials_for_line_3;
    // we put as a precondition that the polynomial for spheres represents
    // a sphere and not an isolated point or an empty_space
    CGAL_kernel_precondition(!(e3.empty_space() || e3.isolated_point())); 
    // if the planes are the same
    // since the solution can only be points or nothing
    // the only case that remains is a plane tangent to a sphere
    if(same_solutions<RT>(e1,e2)) {
      // they are tangent
      return internal::solve_tangent<AK>(e1,e3,res);
    }
    if(!intersect<AK>(e1,e2)) return res;
    Polynomials_for_line_3 pl = line_from_2_planes<AK>(e1,e2);
    return solve<AK>(pl,e3,res);
  }

  template < class AK, class OutputIterator >
  inline 
  OutputIterator
    solve(const typename AK::Polynomial_for_spheres_2_3& e1,
	  const typename AK::Polynomial_1_3 & e2,
	  const typename AK::Polynomial_1_3 & e3,
	 OutputIterator res )
  {
    return solve<AK>(e2,e3,e1,res);
  }

  template < class AK, class OutputIterator >
  inline 
  OutputIterator
    solve(const typename AK::Polynomial_for_spheres_2_3& e1,
	  const typename AK::Polynomial_for_spheres_2_3& e2,
	  const typename AK::Polynomial_1_3 & e3,
	 OutputIterator res )
  {
    typedef typename AK::Polynomial_1_3 Polynomial_1_3;
    // we put as a precondition that the polynomial for spheres represents
    // a sphere and not an isolated point or an empty_space
    CGAL_kernel_precondition(!(e1.empty_space() || e1.isolated_point())); 
    CGAL_kernel_precondition(!(e2.empty_space() || e2.isolated_point())); 
    // The solve can only be points or nothing
    if(e1 == e2) {
      if(tangent<AK>(e3,e1)) {
        return internal::solve_tangent<AK>(e3,e1,res);
      }
      CGAL_kernel_precondition(!(intersect<AK>(e3,e1)));
      return res;
    }
    if(intersect<AK>(e1,e2)) {
      const Polynomial_1_3 p1 = plane_from_2_spheres<AK>(e1,e2);
      return solve<AK>(p1,e3,e1,res);
    } return res;
  }

  template < class AK, class OutputIterator >
  inline 
  OutputIterator
    solve(const typename AK::Polynomial_1_3 & e1,
          const typename AK::Polynomial_for_spheres_2_3& e2,
	  const typename AK::Polynomial_for_spheres_2_3& e3,	  
	 OutputIterator res )
  {
    return solve<AK>(e2,e3,e1,res);
  }

  template < class AK, class OutputIterator >
  inline 
  OutputIterator
    solve( const std::pair<typename AK::Polynomial_for_spheres_2_3, typename AK::Polynomial_1_3 > & e1,
	 const typename AK::Polynomial_1_3 & e2,
	 OutputIterator res )
  {
    return solve<AK>(e1.first,e1.second,e2,res);
  }

  template < class AK, class OutputIterator >
  inline 
  OutputIterator
    solve(const typename AK::Polynomial_1_3 & e1, 
          const std::pair<typename AK::Polynomial_for_spheres_2_3, typename AK::Polynomial_1_3 > & e2,
	  OutputIterator res )
  {
    return solve<AK>(e2,e1,res);
  }

  template < class AK, class OutputIterator >
  inline 
  OutputIterator
    solve( const std::pair<typename AK::Polynomial_for_spheres_2_3, typename AK::Polynomial_1_3 > & e1,
	 const typename AK::Polynomial_for_spheres_2_3 & e2,
	 OutputIterator res )
  {
    return solve<AK>(e1.first, e2, e1.second,res);
  }

  template < class AK, class OutputIterator >
  inline 
  OutputIterator
    solve(const typename AK::Polynomial_for_spheres_2_3 & e1, 
          const std::pair<typename AK::Polynomial_for_spheres_2_3, typename AK::Polynomial_1_3 > & e2,
	  OutputIterator res )
  {
    return solve<AK>(e2,e1, res);
  }

  template < class AK, class OutputIterator >
  inline 
  OutputIterator
    solve( const std::pair<typename AK::Polynomial_for_spheres_2_3, typename AK::Polynomial_1_3 > & e1,
	 const std::pair<typename AK::Polynomial_for_spheres_2_3, typename AK::Polynomial_1_3 > & e2,
	 OutputIterator res )
  {
    typedef typename AK::RT RT;
    typedef typename AK::Root_for_spheres_2_3 Root_for_spheres_2_3;
    typedef typename AK::Polynomial_for_spheres_2_3 Polynomial_for_spheres_2_3;
    typedef typename AK::Polynomial_1_3 Polynomial_1_3;
    const Polynomial_for_spheres_2_3 &s1 = e1.first;
    const Polynomial_for_spheres_2_3 &s2 = e2.first;
    const Polynomial_1_3 &p1 = e1.second;
    const Polynomial_1_3 &p2 = e2.second;

    // The two circles cannot be equal
    // this equality function assumes that the polynomial for sphere represent
    // the sphere with the least radius that produces the circle
    // this is an unique representation, and it is more efficient
    // since the function to see if 2 circles are equal are much more efficient
    CGAL_kernel_precondition((!(e1 == e2)) || (same_solutions<RT>(p1,p2)));

    // Let say that the 2 pairs of equation must be a circle or a point at least
    // ASK: What about removing this pre-condition since it is an Algebraic Kernel?
    CGAL_kernel_precondition((intersect<AK>(p1,s1)));
    CGAL_kernel_precondition((intersect<AK>(p2,s2)));

    if(p1.empty_space()) return res;
    if(p2.empty_space()) return res;
    if(p1.undefined()) {
      return solve<AK>(s1, s2, p2, res);
    }
    if(p2.undefined()) {
      return solve<AK>(s1, s2, p1, res);
    }
    if(same_solutions<RT>(p1, p2)) { 
      return solve<AK>(s1, s2, p1, res);
    }

    typedef std::vector< std::pair<Root_for_spheres_2_3, int> > solutions_container;
    solutions_container solutions;
    solve<AK>(p1, p2, s1, std::back_inserter(solutions));
    if(solutions.size() == 0) return res;
    if(solutions.size() == 1) {
      if(sign_at<AK>(s2, solutions[0].first) == ZERO) {
        *res++ = solutions[0]; 
      } return res;
    }

    // number of solution = 2, we need to set the correct multiplicity
    bool k1 = (sign_at<AK>(s2, solutions[0].first) == ZERO),
         k2 = (sign_at<AK>(s2, solutions[1].first) == ZERO);
    if(k1 && k2) {
      *res++ = solutions[0];
      *res++ = solutions[1];
      return res;
    }
    if(k1) {
      solutions[0].second = 2u;
      *res++ = solutions[0];
      return res;
    }
    if(k2) {
      solutions[1].second = 2u;
      *res++ = solutions[1];
      return res;
    }
    return res;
  }

template < class AK, class OutputIterator >
  inline 
  OutputIterator
    solve( const std::pair<typename AK::Polynomial_for_spheres_2_3, typename AK::Polynomial_1_3 > & e1,
	   const typename AK::Polynomials_for_line_3 & l,
	   OutputIterator res )
  {
    typedef typename AK::Root_for_spheres_2_3 Root_for_spheres_2_3;
    typedef typename AK::Polynomial_for_spheres_2_3 Polynomial_for_spheres_2_3;
    typedef typename AK::Polynomial_1_3 Polynomial_1_3;
    const Polynomial_for_spheres_2_3 &s1 = e1.first;
    const Polynomial_1_3 &p1 = e1.second;

    // Let say that the 2 pairs of equation must be a circle or a point at least
    // ASK: What about removing this pre-condition since it is an Algebraic Kernel?
    CGAL_kernel_precondition((intersect<AK>(p1,s1)));
    CGAL_kernel_precondition(!(l.degenerated()));

    if(p1.empty_space()) return res;
    if(p1.undefined()) {
      return solve<AK>(s1, l, res);
    }

    typedef std::vector< std::pair<Root_for_spheres_2_3, int> > solutions_container;
    solutions_container solutions;
    solve<AK>(s1, l, std::back_inserter(solutions));
    if(solutions.size() == 0) return res;
    if(solutions.size() == 1) {
      if(sign_at<AK>(p1, solutions[0].first) == ZERO) {
        *res++ = solutions[0]; 
      } return res;
    }

    // number of solution = 2, we need to set the correct multiplicity
    bool k1 = (sign_at<AK>(p1, solutions[0].first) == ZERO),
         k2 = (sign_at<AK>(p1, solutions[1].first) == ZERO);
    if(k1 && k2) {
      *res++ = solutions[0];
      *res++ = solutions[1];
      return res;
    }
    if(k1) {
      solutions[0].second = 2u;
      *res++ = solutions[0];
      return res;
    }
    if(k2) {
      solutions[1].second = 2u;
      *res++ = solutions[1];
      return res;
    }
    return res;
  }

  template < class AK, class OutputIterator >
  inline 
  OutputIterator
    solve( const typename AK::Polynomials_for_line_3 & l,
	   const std::pair<typename AK::Polynomial_for_spheres_2_3, typename AK::Polynomial_1_3 > & e1,
	   OutputIterator res ) {
    return solve<AK>(e1, l, res);
  }

  } // namespace AlgebraicSphereFunctors
} // namespace CGAL

#endif //  CGAL_ALGEBRAIC_KERNEL_FUNCTIONS_ON_ROOTS_AND_POLYNOMIAL_1_3_AND_2_3_H
