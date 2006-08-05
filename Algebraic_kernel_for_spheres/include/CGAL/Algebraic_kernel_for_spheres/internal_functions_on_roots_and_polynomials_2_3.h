// Copyright (c) 2005  INRIA Sophia-Antipolis (France) 
// All rights reserved.
//
// Authors : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//           Sylvain Pion     <Sylvain.Pion@sophia.inria.fr>
//           Julien Hazebrouck
//           Damien Leroy
// 
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

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
        return CGALi::solve_tangent<AK>(p,e1,res);
      }
      CGAL_kernel_precondition(!(intersect<AK>(e1,e3)));
      return res;
    }
    if((e1 == e3) || (e2 == e3)) {
      if(tangent<AK>(e1,e2)) {
        Polynomial_1_3 p = plane_from_2_spheres<AK>(e1,e2);
        return CGALi::solve_tangent<AK>(p,e1,res);
      }
      CGAL_kernel_precondition(!(intersect<AK>(e1,e2)));
      return res;
    }
    
    // non degenerated case
    if(intersect<AK>(e1,e2)) {
      Polynomial_1_3 p1 = plane_from_2_spheres<AK>(e1,e2);
      if(intersect<AK>(e2,e3)) {
        Polynomial_1_3 p2 = plane_from_2_spheres<AK>(e2,e3); 
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

    return Root_for_spheres_2_3(make_root_of_2(s.a(),i?-1:1,s.r_sq()),
                                Root_of_2(s.b()),
                                Root_of_2(s.c()));
  }

  template <class AK, class OutputIterator>
  OutputIterator
  x_critical_points(const typename AK::Polynomial_for_spheres_2_3 & s, OutputIterator res)
  {
    typedef typename AK::Root_of_2 Root_of_2;
    typedef typename AK::Root_for_spheres_2_3 Root_for_spheres_2_3; 

    *res++ =  Root_for_spheres_2_3(make_root_of_2(s.a(),-1,s.r_sq()),
                                Root_of_2(s.b()),
                                Root_of_2(s.c()));
    *res++ =  Root_for_spheres_2_3(make_root_of_2(s.a(),1,s.r_sq()),
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

    return Root_for_spheres_2_3(Root_of_2(s.a()),
                                make_root_of_2(s.b(),i?-1:1,s.r_sq()),
                                Root_of_2(s.c()));
  }
  
  template <class AK, class OutputIterator>
  OutputIterator
  y_critical_points(const typename AK::Polynomial_for_spheres_2_3 & s, OutputIterator res)
  {
    typedef typename AK::Root_of_2 Root_of_2;
    typedef typename AK::Root_for_spheres_2_3 Root_for_spheres_2_3; 

    *res++ =  Root_for_spheres_2_3(Root_of_2(s.a()),
                                make_root_of_2(s.b(),-1,s.r_sq()),
                                Root_of_2(s.c()));
    *res++ =  Root_for_spheres_2_3(Root_of_2(s.a()),
                                make_root_of_2(s.b(),1,s.r_sq()),
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
                                make_root_of_2(s.c(),i?-1:1,s.r_sq()));
  }
  
  template <class AK, class OutputIterator>
  OutputIterator
  z_critical_points(const typename AK::Polynomial_for_spheres_2_3 & s, OutputIterator res)
  {
    typedef typename AK::Root_of_2 Root_of_2;
    typedef typename AK::Root_for_spheres_2_3 Root_for_spheres_2_3; 

    *res++ =  Root_for_spheres_2_3(Root_of_2(s.a()),
                                Root_of_2(s.b()),
                                make_root_of_2(s.c(),-1,s.r_sq()));
    *res++ =  Root_for_spheres_2_3(Root_of_2(s.a()),
                                Root_of_2(s.b()),
                                make_root_of_2(s.c(),1,s.r_sq()));
    return res;
  }

  // ONLY the x is given (Root_of_2), because the other coordinates may be Root_of_4
  template <class AK>
  typename AK::Root_of_2
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
    CGAL_kernel_precondition(sign(p.a() * s.a() + p.b() * s.b() + 
                                  p.c() * s.c() + p.d()) == ZERO);

    const FT sqbc = CGAL::square(p.b()) + CGAL::square(p.c());
    const FT sq_sum = sqbc + CGAL::square(p.a());
    const FT delta = (sqbc * s.r_sq())/sq_sum;

    return make_root_of_2(s.a(),i?-1:1,delta);
  }

  // ONLY the x is given (Root_of_2), because the other coordinates may be Root_of_4
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
    CGAL_kernel_precondition(sign(p.a() * s.a() + p.b() * s.b() + 
                                  p.c() * s.c() + p.d()) == ZERO);

    const FT sqbc = CGAL::square(p.b()) + CGAL::square(p.c());
    const FT sq_sum = sqbc + CGAL::square(p.a());
    const FT delta = (sqbc * s.r_sq())/sq_sum;
    
    *res++ =  make_root_of_2(s.a(),-1,delta);
    *res++ =  make_root_of_2(s.a(),1,delta);
    return res;
  }

  // ONLY the y is given (Root_of_2), because the other coordinates may be Root_of_4
  template <class AK>
  typename AK::Root_of_2
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
    CGAL_kernel_precondition(sign(p.a() * s.a() + p.b() * s.b() + 
                                  p.c() * s.c() + p.d()) == ZERO);

    const FT sqac = CGAL::square(p.a()) + CGAL::square(p.c());
    const FT sq_sum = sqac + CGAL::square(p.b());
    const FT delta = (sqac * s.r_sq())/sq_sum;

    return make_root_of_2(s.b(),i?-1:1,delta);
  }

  // ONLY the y is given (Root_of_2), because the other coordinates may be Root_of_4
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
    CGAL_kernel_precondition(sign(p.a() * s.a() + p.b() * s.b() + 
                                  p.c() * s.c() + p.d()) == ZERO);

    const FT sqac = CGAL::square(p.a()) + CGAL::square(p.c());
    const FT sq_sum = sqac + CGAL::square(p.b());
    const FT delta = (sqac * s.r_sq())/sq_sum;
    
    *res++ =  make_root_of_2(s.b(),-1,delta);
    *res++ =  make_root_of_2(s.b(),1,delta);
    return res;
  }

  // ONLY the z is given (Root_of_2), because the other coordinates may be Root_of_4
  template <class AK>
  typename AK::Root_of_2
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
    CGAL_kernel_precondition(sign(p.a() * s.a() + p.b() * s.b() + 
                                  p.c() * s.c() + p.d()) == ZERO);

    const FT sqab = CGAL::square(p.a()) + CGAL::square(p.b());
    const FT sq_sum = sqab + CGAL::square(p.c());
    const FT delta = (sqab * s.r_sq())/sq_sum;

    return make_root_of_2(s.c(),i?-1:1,delta);
  }

  // ONLY the z is given (Root_of_2), because the other coordinates may be Root_of_4
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
    CGAL_kernel_precondition(sign(p.a() * s.a() + p.b() * s.b() + 
                                  p.c() * s.c() + p.d()) == ZERO);

    const FT sqab = CGAL::square(p.a()) + CGAL::square(p.b());
    const FT sq_sum = sqab + CGAL::square(p.c());
    const FT delta = (sqab * s.r_sq())/sq_sum;
    
    *res++ =  make_root_of_2(s.c(),-1,delta);
    *res++ =  make_root_of_2(s.c(),1,delta);
    return res;
  }

  } // namespace AlgebraicSphereFunctors
} // namespace CGAL

#endif //  CGAL_ALGEBRAIC_KERNEL_FUNCTIONS_ON_ROOTS_AND_POLYNOMIALS_2_3_H
