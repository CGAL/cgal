// ============================================================================
//
// Copyright (c) 2001-2006 Max-Planck-Institut Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of EXACUS (http://www.mpi-inf.mpg.de/projects/EXACUS/).
//
// ----------------------------------------------------------------------------
//
// Library       : CGAL
// File          : test/modular_gcd_utcf_algorithm_M.C
// CGAL_release   : $Name:  $
// Revision      : $Revision$
// Revision_date : $Date$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Dominik Hülse <dominik.huelse@gmx.de>
//                 Michael Hemmer <hemmer@mpi-inf.mpg.de>
//
// ============================================================================
/*! \file CGAL/polynomial_gcd.C
    test for the modular algorithm modular_gcd_utcf_algorithm_M to compute the gcd of
    univariate polynomials with integer coefficients
*/

#define MY_FUNCTION_CALL modular_gcd_utcf_algorithm_M


#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/CORE_arithmetic_kernel.h>
#include <CGAL/LEDA_arithmetic_kernel.h>
#include <CGAL/gen_polynomials.h>
#include <CGAL/test_modular_gcd.h>

int main(){

  // Set wrong rounding mode to test modular arithmetic
  CGAL::Protect_FPU_rounding<true> pfr(CGAL_FE_UPWARD);
  #ifdef CGAL_USE_LEDA
   CGAL::internal::test_modular_gcd<CGAL::LEDA_arithmetic_kernel>
       (CGAL::Unique_factorization_domain_tag());
  #endif // CGAL_USE_LEDA

  #ifdef CGAL_USE_CORE
   CGAL::internal::test_modular_gcd<CGAL::CORE_arithmetic_kernel>
       (CGAL::Unique_factorization_domain_tag());
  #endif // Lis_HAVE_CORE

  return 0;
}


// EOF
