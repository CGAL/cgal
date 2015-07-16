// ============================================================================
//
// Copyright (c) 2001-2006 Max-Planck-Institut Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of EXACUS (http://www.mpi-inf.mpg.de/projects/EXACUS/).
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
// ----------------------------------------------------------------------------
//
// Library       : CGAL
// File          : test/modular_gcd_utils.C
// CGAL_release   : $Name:  $
// Revision      : $Revision$
// Revision_date : $Date$
//
// Author(s)     : Dominik Huelse <dominik.huelse@gmx.de>
//
//
// ============================================================================

/*! \file CGAL/Polynomial/modular_gcd_utils.C 
  test for the function euclidean_division_obstinate
*/

#include <CGAL/basic.h>
#include <CGAL/gen_polynomials.h>
#include <CGAL/Polynomial/modular_gcd_utils.h>
#include <CGAL/Random.h>
#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/CORE_arithmetic_kernel.h>
#include <CGAL/LEDA_arithmetic_kernel.h>
#include <CGAL/Sqrt_extension.h>
#include <CGAL/Polynomial.h>
#include <CGAL/use.h>

#include <cassert>

template<class AT>
void test_modular_gcd_utils() {
  CGAL_SNAP_ARITHMETIC_KERNEL_TYPEDEFS(AT);
  CGAL_USE_TYPE(Rational);
  CGAL_USE_TYPE(Field_with_sqrt);

  CGAL::Random my_random(4711);
  ::CGAL::set_pretty_mode(std::cout);

  typedef typename AT::Integer Integer;
  typedef Integer NT;

  {  
    // testing integers
    Integer f, g, q, r;
        
    // random integers 
    for(int l=0;l<100;l++){
      f = CGAL::internal::rand_int<Integer>(my_random.get_int(10,1000));
      g = CGAL::internal::rand_int<Integer>(my_random.get_int(10,1000));
      CGAL::internal::euclidean_division_obstinate(f, g, q, r);
      assert(f==g*q+r);
    }
  }
    
  {  
    typedef CGAL::Polynomial<Integer> Poly;
    typedef typename CGAL::Modular_traits<Poly>::Residue_type   MPoly;
    CGAL::Residue::set_current_prime(43);
        
    Poly a, b;
    MPoly mQ, mR;

    a = Poly(0);
    b = Poly(NT(122),NT(72));
  
    MPoly ma = CGAL::modular_image(a);
    MPoly mb = CGAL::modular_image(b);
    CGAL::internal::euclidean_division_obstinate(ma, mb, mQ, mR);
    assert(ma==mb*mQ+mR);
        
    // random polynomials with integer coefficients
    Poly   f, g;
    int l;

    // random polynomials with integer coefficients
    for(l=0;l<100;l++){
      f = CGAL::internal::rand_Poly_int<Integer>(my_random.get_int(10,1000),
          my_random.get_int(1,15));
      g = CGAL::internal::rand_Poly_int<Integer>(my_random.get_int(10,1000), 
          my_random.get_int(1,10));
      CGAL::Residue::set_current_prime(4483);
      MPoly mf= CGAL::modular_image(f);
      MPoly mg= CGAL::modular_image(g);
      CGAL::internal::euclidean_division_obstinate(mf, mg, mQ, mR);

      assert(mf ==mg*mQ+mR);
    }
  }

  { 
    typedef CGAL::Sqrt_extension<Integer,Integer> EXT;
    typedef CGAL::Polynomial<EXT>  Poly_Ext;
    typedef typename CGAL::Modular_traits<Poly_Ext>::Residue_type   MPoly_Ext;
    CGAL::Residue::set_current_prime(43);
        
    Poly_Ext a, b;
    MPoly_Ext mQ, mR;

    a = Poly_Ext(0);
    b = Poly_Ext(NT(122),NT(72),NT(17));
  
    MPoly_Ext ma = CGAL::modular_image(a);
    MPoly_Ext mb = CGAL::modular_image(b);
    CGAL::internal::euclidean_division_obstinate(ma, mb, mQ, mR);
    assert(ma == mb*mQ+mR);
  
    // random polynomials with sqrt coefficients
    Poly_Ext   f, g;
    int l;
    MPoly_Ext mf, mg;

    CGAL::Residue::set_current_prime(4253);
    for(l=0;l<2;l++){
      f = CGAL::internal::rand_Poly_sqrt<EXT,Integer>
        (my_random.get_int(10,1000),my_random.get_int(1,15),NT(8293));
      g = CGAL::internal::rand_Poly_sqrt<EXT,Integer>
        (my_random.get_int(10,1000),my_random.get_int(1,10),NT(8293));
      MPoly_Ext mf= CGAL::modular_image(f);
      MPoly_Ext mg= CGAL::modular_image(g);
      CGAL::internal::euclidean_division_obstinate(mf, mg, mQ, mR);
      assert(mf == mg*mQ+mR);
    }
  }
}
int main(){
  // Enforce IEEE double precision and to nearest before using modular arithmetic
  CGAL::Protect_FPU_rounding<true> pfr(CGAL_FE_TONEAREST);  
 
#ifdef CGAL_USE_LEDA
  test_modular_gcd_utils<CGAL::LEDA_arithmetic_kernel>();       
#endif // CGAL_USE_LEDA    

#ifdef CGAL_USE_CORE      
  test_modular_gcd_utils<CGAL::CORE_arithmetic_kernel>();
#endif // Lis_HAVE_CORE
     
  return 0;
}


// EOF
 
