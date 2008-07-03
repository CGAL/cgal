// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     : Michael Kerber <mkerber@mpi-inf.mpg.de>
//
// ============================================================================

// This file defines several instances of Algebraic_kernel_1 that are often
// used in tests and demos. 

#ifndef CGAL_PREFERED_ALGEBRAIC_KERNEL_H
#define CGAL_PREFERED_ALGEBRAIC_KERNEL_H 1

#include <CGAL/basic.h>

#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Algebraic_kernel_1.h>

// Needed for the "fast" kernel
#include <CGAL/Algebraic_kernel_d/Algebraic_real_quadratic_refinement_rep_bfi.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_descartes.h>

// Needed for the bisection kernel
#include <CGAL/Algebraic_kernel_d/Algebraic_real_pure.h>
#include <CGAL/Algebraic_kernel_d/Descartes.h>

CGAL_BEGIN_NAMESPACE

/**
 * Defines default and a fast Algebraic_kernel_1, 
 * depending on the coefficient type.
 */
template<typename Coefficient, 
         typename Boundary = typename CGAL::Get_arithmetic_kernel< Coefficient >::Arithmetic_kernel::Rational>
struct Get_algebraic_kernel_1 {

    typedef CGAL::Algebraic_kernel_1 <Coefficient, Boundary>
        Default_algebraic_kernel_1;

    typedef CGAL::Algebraic_kernel_1
    < Coefficient, 
      Boundary,
      CGAL::CGALi::Algebraic_real_rep< Coefficient, Boundary >,
      CGAL::CGALi::Descartes< CGAL::Polynomial< Coefficient >, Boundary >
    > Algebraic_kernel_with_bisection_and_descartes_1;

    typedef CGAL::Algebraic_kernel_1
    < Coefficient, 
      Boundary,
      CGAL::CGALi::Algebraic_real_quadratic_refinement_rep_bfi
           < Coefficient, Boundary >,
      CGAL::CGALi::Bitstream_descartes
        < CGAL::Polynomial< Coefficient >, Boundary >
    > Algebraic_kernel_with_qir_and_bitstream_1;
};

CGAL_END_NAMESPACE

#endif
