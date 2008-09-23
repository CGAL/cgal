// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     : Sebastian Limbach <slimbach@mpi-inf.mpg.de>
//                 Michael Hemmer    <hemmer@mpi-inf.mpg.de>
//
// ============================================================================

// Test of Algebraic_kernel

#include <CGAL/basic.h>
#include <CGAL/Algebraic_kernel_1.h>
#include <CGAL/Algebraic_kernel_d/Algebraic_real_rep_bfi.h>
#include <CGAL/Algebraic_kernel_d/Algebraic_real_rep.h>
#include <CGAL/Algebraic_kernel_d/Algebraic_real_quadratic_refinement_rep_bfi.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_descartes.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_descartes_rndl_tree_traits.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_coefficient_kernel.h>
#include <CGAL/Algebraic_kernel_d/Descartes.h>
#include <CGAL/_test_algebraic_kernel_1.h>

#include <CGAL/Arithmetic_kernel.h>


template< class Coefficient_, class Boundary_, class RepClass >
void test_algebraic_kernel_coeff_boundary_rep() {
  typedef Coefficient_ Coefficient;
  typedef Boundary_    Boundary;
  typedef RepClass     Rep_class;
  
  typedef typename CGAL::Polynomial_type_generator<Coefficient,1>::Type 
      Polynomial_1;
  typedef CGAL::CGALi::Algebraic_real_pure
    < Coefficient, Boundary, CGAL::Handle_policy_no_union, Rep_class >   Algebraic_real_1;
  
  typedef CGAL::CGALi::Descartes< Polynomial_1, Boundary >               Descartes;
  typedef CGAL::CGALi::Bitstream_descartes<
            CGAL::CGALi::Bitstream_descartes_rndl_tree_traits
            <CGAL::CGALi::Bitstream_coefficient_kernel<Coefficient> > >  BDescartes;
//   typedef CGAL::CGALi::Bitstream_descartes< Polynomial_1, Boundary > BDescartes;
  
  typedef CGAL::Algebraic_kernel_1< Coefficient, Boundary, Rep_class , Descartes>
    Kernel_Descartes;
  typedef CGAL::Algebraic_kernel_1< Coefficient, Boundary, Rep_class , BDescartes>
    Kernel_BDescartes;
  
  CGAL::CGALi::test_algebraic_kernel_1< 
    Kernel_Descartes, Algebraic_real_1, Descartes, Coefficient, Polynomial_1, Boundary >();   
  CGAL::CGALi::test_algebraic_kernel_1< 
    Kernel_BDescartes, Algebraic_real_1, BDescartes, Coefficient, Polynomial_1, Boundary >();    

}


template< class Coeff, class Boundary >
void test_algebraic_kernel_coeff_boundary() {
  test_algebraic_kernel_coeff_boundary_rep<Coeff,Boundary,
    CGAL::CGALi::Algebraic_real_rep< Coeff, Boundary > > ();
  test_algebraic_kernel_coeff_boundary_rep<Coeff,Boundary,
    CGAL::CGALi::Algebraic_real_rep_bfi< Coeff, Boundary > > ();
  test_algebraic_kernel_coeff_boundary_rep<Coeff,Boundary,
    CGAL::CGALi::Algebraic_real_quadratic_refinement_rep_bfi< Coeff, Boundary > > ();
}


template< class ArithmeticKernel >
void test_algebraic_kernel() {
    typedef ArithmeticKernel AK;
    typedef typename AK::Integer Integer;
    typedef typename AK::Rational Rational; 

    test_algebraic_kernel_coeff_boundary<Integer, Rational>();
    test_algebraic_kernel_coeff_boundary<Rational, Rational>();
    test_algebraic_kernel_coeff_boundary
      <CGAL::Sqrt_extension< Integer, Integer>, Rational>();
    test_algebraic_kernel_coeff_boundary
      <CGAL::Sqrt_extension< Rational, Integer>, Rational>();
    test_algebraic_kernel_coeff_boundary
      <CGAL::Sqrt_extension< Rational, Rational>, Rational>();
}

int main() {
#ifdef CGAL_USE_LEDA
    test_algebraic_kernel< CGAL::LEDA_arithmetic_kernel >();
#endif
#ifdef CGAL_USE_CORE
    test_algebraic_kernel< CGAL::CORE_arithmetic_kernel >();
#endif
    
    return 0;
}
