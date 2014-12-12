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
// File          : test/modular_gcd_utcf_dfai.C
// CGAL_release   : $Name:  $
// Revision      : $Revision$
// Revision_date : $Date$
//
// Author(s)     : Dominik Hülse <dominik.huelse@gmx.de>
//                 Michael Hemmer <hemmer@mpi-inf.mpg.de>
//
// ============================================================================

/*! \file CGAL/polynomial_gcd.C 
    test for the modular algorithm modular_gcd_utcf_dfai to compute the gcd of 
    univariate polynomials with (nested) sqrt-extensions
*/

#define MY_FUNCTION_CALL modular_gcd_utcf_dfai

#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/CORE_arithmetic_kernel.h>
#include <CGAL/LEDA_arithmetic_kernel.h>
#include <CGAL/test_modular_gcd.h>

int main(){
    
  // Set wrong rounding mode to test modular arithmetic 
  CGAL::Protect_FPU_rounding<true> pfr(CGAL_FE_UPWARD);

  #ifdef CGAL_USE_LEDA
   CGAL::internal::test_modular_gcd<CGAL::LEDA_arithmetic_kernel>
       (CGAL::Integral_domain_tag());
  #endif // CGAL_USE_LEDA    

  #ifdef CGAL_USE_CORE      
   CGAL::internal::test_modular_gcd<CGAL::CORE_arithmetic_kernel>
       (CGAL::Integral_domain_tag());
  #endif // Lis_HAVE_CORE
     
  return 0;
}


// EOF
 
