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
//             Julien Hazebrouck
//             Damien Leroy

#ifndef CGAL_ROOT_FOR_SPHERES_2_3_H
#define CGAL_ROOT_FOR_SPHERES_2_3_H

#include <CGAL/license/Circular_kernel_3.h>


#include <iostream>
#include <CGAL/Polynomials_1_3.h>
#include <CGAL/Polynomials_2_3.h>
#include <CGAL/Polynomials_for_line_3.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Root_of_traits.h>

namespace CGAL {

template < typename RT_ >
class Root_for_spheres_2_3 {

  typedef RT_                                        RT;
  typedef typename Root_of_traits< RT >::RootOf_2    Root_of_2;
  typedef typename Root_of_traits< RT >::RootOf_1    FT;
  typedef CGAL::Polynomial_1_3< FT >                 Polynomial_1_3;
  typedef CGAL::Polynomial_for_spheres_2_3< FT >     Polynomial_for_spheres_2_3;
  typedef CGAL::Polynomials_for_line_3< FT >         Polynomials_for_line_3;

  private:
    Root_of_2 x_;
    Root_of_2 y_;
    Root_of_2 z_;
    
  public:
  Root_for_spheres_2_3(){}
  
  
  Root_for_spheres_2_3(const Root_of_2& r1,
		       const Root_of_2& r2,
		       const Root_of_2& r3)
    : x_(r1), y_(r2), z_(r3)
  {
    // This assertion sont work if Root_of_2 is 
    // Interval_nt (and dont have is_rational, gamma, etc..)
    /*CGAL_assertion(
                ((r1.is_rational() && r2.is_rational()) ||
                 (r1.is_rational() && r3.is_rational()) ||
                 (r2.is_rational() && r3.is_rational()) ||
                 ((r1.is_rational()) && (r2.gamma() == r3.gamma())) ||
                 ((r2.is_rational()) && (r1.gamma() == r3.gamma())) ||
                 ((r3.is_rational()) && (r1.gamma() == r2.gamma())) ||
                 ((r1.gamma() == r2.gamma()) && (r2.gamma() == r3.gamma())))
    );*/
  }

  const Root_of_2& x() const 
  { return x_; }
    
  const Root_of_2& y() const 
  { return y_; }

  const Root_of_2& z() const 
  { return z_; }

  // On fait l'evaluation de (x,y,z) pour le plan 
  // aX + bY + cZ + d, donne
  const Root_of_2 evaluate(const Polynomial_1_3 &p) const {
    return (p.a() * x()) + (p.b() * y()) + (p.c() * z()) + p.d();
  }

  // On fait l'evaluation de (x,y,z) pour le plan 
  // (X-a)^2 + (Y-b)^2 + (Z-c)^2 - r_sq, donne
  const Root_of_2 evaluate(const Polynomial_for_spheres_2_3 &p) const {
    return square(x() - p.a()) +
           square(y() - p.b()) +
           square(z() - p.c()) -
           p.r_sq();
  }

  // On verifie si (x,y,z) fait partie la ligne donne
  bool is_on_line(const Polynomials_for_line_3 &p) const {
    Root_of_2 t;
    bool already = false;
    if(!is_zero(p.a1())) {
      t = (x() - p.b1())/p.a1(); 
      already = true;
    } else if(p.b1() != x()) return false;
    if(!is_zero(p.a2())) {
      if(!already) {
        t = (y() - p.b2())/p.a2(); 
        already = true;
      }
      else if((p.a2() * t + p.b2()) != y()) return false;
    } else if(p.b2() != y()) return false;
    if(!is_zero(p.a3())) {
      if(!already) return true;
      else if((p.a3() * t + p.b3()) != z()) return false;
    } else if(p.b3() != z()) return false;
    return true;
  }

  CGAL::Bbox_3 bbox() const
  {
    const Root_of_2 &ox = x();
    const Root_of_2 &oy = y();
    const Root_of_2 &oz = z();

    CGAL::Interval_nt<> 
        ix=to_interval(ox),
        iy=to_interval(oy),
        iz=to_interval(oz);
      return CGAL::Bbox_3(ix.inf(),iy.inf(),iz.inf(),
	                ix.sup(),iy.sup(),iz.sup());
    /* 
    // Note: This is a more efficient version
    // but it won't work (in the future) 
    // with some Lazy_Curved_kernel_3
    // because is_rational(), gamma(), etc.. is not defined
    // for Interval_nt<false> data type	
    const Root_of_2 &ox = x();
    const Root_of_2 &oy = y();
    const Root_of_2 &oz = z();

    const bool x_rat = ox.is_rational();
    const bool y_rat = oy.is_rational();
    const bool z_rat = oz.is_rational();

    if(((x_rat?1:0) + (y_rat?1:0) +(z_rat?1:0)) > 1) {
      CGAL::Interval_nt<> 
        ix=to_interval(ox),
        iy=to_interval(oy),
        iz=to_interval(oz);
      return CGAL::Bbox_3(ix.inf(),iy.inf(),iz.inf(),
	                ix.sup(),iy.sup(),iz.sup());
    }

    if(z_rat) {
      const CGAL::Interval_nt<true> alpha1 = to_interval(ox.alpha());
      const CGAL::Interval_nt<true> beta1 = to_interval(ox.beta());
      const CGAL::Interval_nt<true> alpha2 = to_interval(oy.alpha());
      const CGAL::Interval_nt<true> beta2 = to_interval(oy.beta());
      const CGAL::Interval_nt<true> g = to_interval(ox.gamma());
      const CGAL::Interval_nt<true> sqrtg = CGAL::sqrt(g);
      const CGAL::Interval_nt<true> ix = alpha1 + beta1 * sqrtg;
      const CGAL::Interval_nt<true> iy = alpha2 + beta2 * sqrtg;
      const CGAL::Interval_nt<true> iz = to_interval(oz);
      return CGAL::Bbox_3(ix.inf(),iy.inf(),iz.inf(),
	                ix.sup(),iy.sup(),iz.sup());
    }

    if(y_rat) {
      const CGAL::Interval_nt<true> alpha1 = to_interval(ox.alpha());
      const CGAL::Interval_nt<true> beta1 = to_interval(ox.beta());
      const CGAL::Interval_nt<true> alpha2 = to_interval(oz.alpha());
      const CGAL::Interval_nt<true> beta2 = to_interval(oz.beta());
      const CGAL::Interval_nt<true> g = to_interval(ox.gamma());
      const CGAL::Interval_nt<true> sqrtg = CGAL::sqrt(g);
      const CGAL::Interval_nt<true> ix = alpha1 + beta1 * sqrtg;
      const CGAL::Interval_nt<true> iz = alpha2 + beta2 * sqrtg;
      const CGAL::Interval_nt<true> iy = to_interval(oy);
      return CGAL::Bbox_3(ix.inf(),iy.inf(),iz.inf(),
	                ix.sup(),iy.sup(),iz.sup());
    }

    if(x_rat) {
      const CGAL::Interval_nt<true> alpha1 = to_interval(oy.alpha());
      const CGAL::Interval_nt<true> beta1 = to_interval(oy.beta());
      const CGAL::Interval_nt<true> alpha2 = to_interval(oz.alpha());
      const CGAL::Interval_nt<true> beta2 = to_interval(oz.beta());
      const CGAL::Interval_nt<true> g = to_interval(oy.gamma());
      const CGAL::Interval_nt<true> sqrtg = CGAL::sqrt(g);
      const CGAL::Interval_nt<true> iy = alpha1 + beta1 * sqrtg;
      const CGAL::Interval_nt<true> iz = alpha2 + beta2 * sqrtg;
      const CGAL::Interval_nt<true> ix = to_interval(ox);
      return CGAL::Bbox_3(ix.inf(),iy.inf(),iz.inf(),
	                ix.sup(),iy.sup(),iz.sup());
    }

    const CGAL::Interval_nt<true> alpha1 = to_interval(ox.alpha());
    const CGAL::Interval_nt<true> beta1 = to_interval(ox.beta());
    const CGAL::Interval_nt<true> alpha2 = to_interval(oy.alpha());
    const CGAL::Interval_nt<true> beta2 = to_interval(oy.beta());
    const CGAL::Interval_nt<true> alpha3 = to_interval(oz.alpha());
    const CGAL::Interval_nt<true> beta3 = to_interval(oz.beta());
    const CGAL::Interval_nt<true> g = to_interval(ox.gamma());
    const CGAL::Interval_nt<true> sqrtg = CGAL::sqrt(g);
    const CGAL::Interval_nt<true> ix = alpha1 + beta1 * sqrtg;
    const CGAL::Interval_nt<true> iy = alpha2 + beta2 * sqrtg;
    const CGAL::Interval_nt<true> iz = alpha3 + beta3 * sqrtg; 
    return CGAL::Bbox_3(ix.inf(),iy.inf(),iz.inf(),
	                ix.sup(),iy.sup(),iz.sup());
    */
  }

};

template < typename RT >
bool 
operator == ( const Root_for_spheres_2_3<RT>& r1,
	      const Root_for_spheres_2_3<RT>& r2 )
{ return (r1.x() == r2.x()) && (r1.y() == r2.y()) && (r1.z() == r2.z()); }

template < typename RT >
std::ostream &
operator<<(std::ostream & os, const Root_for_spheres_2_3<RT> &r)
{ return os << r.x() << " " << r.y() << " " << r.z() << " "; }

template < typename RT >
std::istream &
operator>>(std::istream & is, Root_for_spheres_2_3<RT> &r)
{
  typedef typename Root_of_traits< RT >::RootOf_2         Root_of_2;
  Root_of_2 x,y,z;
  
  is >> x >> y >> z;
  if(is)
    r = Root_for_spheres_2_3<RT>(x,y,z);
  return is;
}

} //namespace CGAL

#endif // CGAL_ROOT_FOR_SPHERES_2_3_H
