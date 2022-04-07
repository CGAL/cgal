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

#define CGAL_TEST_ALL_AK_VARIANTS 1

#include <CGAL/Algebraic_kernel_d_1.h>
#include <CGAL/Algebraic_kernel_d/Algebraic_real_rep_bfi.h>
#include <CGAL/Algebraic_kernel_d/Algebraic_real_rep.h>
#include <CGAL/Algebraic_kernel_d/Algebraic_real_quadratic_refinement_rep_bfi.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_descartes.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_descartes_rndl_tree_traits.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_coefficient_kernel.h>
#include <CGAL/Algebraic_kernel_d/Descartes.h>
#include <CGAL/_test_algebraic_kernel_1.h>

#include <CGAL/Arithmetic_kernel.h>


template< class Coefficient_, class Bound_, class RepClass >
void test_algebraic_kernel_coeff_bound_rep() {
  typedef Coefficient_  Coefficient;
  typedef Bound_        Bound;
  typedef RepClass      Rep_class;

  typedef typename CGAL::Polynomial_type_generator<Coefficient,1>::Type
    Polynomial_1;
  typedef CGAL::internal::Algebraic_real_d_1
    < Coefficient, Bound, CGAL::Handle_policy_no_union, Rep_class >   Algebraic_real_1;
  CGAL_USE_TYPE(Algebraic_real_1);

  typedef CGAL::internal::Descartes< Polynomial_1, Bound >               Descartes;
  typedef CGAL::internal::Bitstream_descartes<
  CGAL::internal::Bitstream_descartes_rndl_tree_traits
    <CGAL::internal::Bitstream_coefficient_kernel<Coefficient> > >  BDescartes;
//   typedef CGAL::internal::Bitstream_descartes< Polynomial_1, Bound > BDescartes;

  typedef CGAL::Algebraic_kernel_d_1< Coefficient, Bound, Rep_class , Descartes>
    Kernel_Descartes;
  typedef CGAL::Algebraic_kernel_d_1< Coefficient, Bound, Rep_class , BDescartes>
    Kernel_BDescartes;

  CGAL::test_algebraic_kernel_1(Kernel_Descartes());
#if CGAL_TEST_ALL_AK_VARIANTS
  CGAL::test_algebraic_kernel_1(Kernel_BDescartes());
#endif
}


template< class Coeff, class Bound >
void test_algebraic_kernel_coeff_bound() {
  test_algebraic_kernel_coeff_bound_rep<Coeff,Bound,
    CGAL::internal::Algebraic_real_rep< Coeff, Bound > > ();
#if CGAL_TEST_ALL_AK_VARIANTS
  test_algebraic_kernel_coeff_bound_rep<Coeff,Bound,
    CGAL::internal::Algebraic_real_rep_bfi< Coeff, Bound > > ();
  test_algebraic_kernel_coeff_bound_rep<Coeff,Bound,
    CGAL::internal::Algebraic_real_quadratic_refinement_rep_bfi< Coeff, Bound > > ();
#endif
}


template< class ArithmeticKernel >
void test_algebraic_kernel() {
  typedef ArithmeticKernel AK;
  typedef typename AK::Integer Integer;
  typedef typename AK::Rational Rational;

  test_algebraic_kernel_coeff_bound<Integer, Rational>();
#if CGAL_TEST_ALL_AK_VARIANTS
  test_algebraic_kernel_coeff_bound<Rational, Rational>();
  test_algebraic_kernel_coeff_bound
    <CGAL::Sqrt_extension< Integer, Integer>, Rational>();
  test_algebraic_kernel_coeff_bound
    <CGAL::Sqrt_extension< Rational, Integer>, Rational>();
  test_algebraic_kernel_coeff_bound
    <CGAL::Sqrt_extension< Rational, Rational>, Rational>();
#endif
}
