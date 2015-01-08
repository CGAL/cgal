#include <iostream>
#include <CGAL/basic.h>
#include <CGAL/use.h>
#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Quotient.h>
#include <CGAL/Test/_test_algebraic_structure.h>
#include <CGAL/Test/_test_real_embeddable.h>
#include <CGAL/Test/_test_fraction_traits.h>
#include <CGAL/Test/_test_rational_traits.h>

template< class AT >
void test_quotient() {
  {
    typedef CGAL::Quotient< typename AT::Integer > NT;
    typedef CGAL::Field_tag Tag;
    typedef CGAL::Tag_true Is_exact;
    CGAL::test_algebraic_structure<NT,Tag, Is_exact>();
    CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(4),NT(6),NT(15));
    CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(-4),NT(6),NT(15));
    CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(4),NT(-6),NT(15));
    CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(-4),NT(-6),NT(15));
    CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(4),NT(6),NT(-15));
    CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(-4),NT(6), NT(15));
    CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(4),NT(-6),NT(-15));
    CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(-4),NT(-6),NT(-15));

    CGAL::test_algebraic_structure<NT,Tag, Is_exact>( NT(5,74), NT(3,25), NT(7,3));
    CGAL::test_algebraic_structure<NT,Tag, Is_exact>(-NT(5,74), NT(3,25), NT(7,3));
    CGAL::test_algebraic_structure<NT,Tag, Is_exact>( NT(5,74),-NT(3,25), NT(7,3));
    CGAL::test_algebraic_structure<NT,Tag, Is_exact>(-NT(5,74),-NT(3,25), NT(7,3));
    CGAL::test_algebraic_structure<NT,Tag, Is_exact>( NT(5,74), NT(3,25),-NT(7,3));
    CGAL::test_algebraic_structure<NT,Tag, Is_exact>(-NT(5,74), NT(3,25),-NT(7,3));
    CGAL::test_algebraic_structure<NT,Tag, Is_exact>( NT(5,74),-NT(3,25),-NT(7,3));
    CGAL::test_algebraic_structure<NT,Tag, Is_exact>(-NT(5,74),-NT(3,25),-NT(7,3));

    CGAL::test_real_embeddable<NT>();
    CGAL::test_fraction_traits<NT>(); 
    // backward compatiblity
    CGAL::test_rational_traits<NT>();  

  }
  /* // Quotient for inexact types not implemented 
     {
      typedef CGAL::Quotient< leda_bigfloat > NT;
      typedef CGAL::Field_with_sqrt_tag Tag;
      CGAL::test_algebraic_structure<NT,Tag>();
      CGAL::test_algebraic_structure<NT,Tag>( NT(5,74), NT(3,25), NT(7,3));
      CGAL::test_algebraic_structure<NT,Tag>(-NT(5,74), NT(3,25), NT(7,3));
      CGAL::test_algebraic_structure<NT,Tag>( NT(5,74),-NT(3,25), NT(7,3));
      CGAL::test_algebraic_structure<NT,Tag>(-NT(5,74),-NT(3,25), NT(7,3));
      CGAL::test_algebraic_structure<NT,Tag>( NT(5,74), NT(3,25),-NT(7,3));
      CGAL::test_algebraic_structure<NT,Tag>(-NT(5,74), NT(3,25),-NT(7,3));
      CGAL::test_algebraic_structure<NT,Tag>( NT(5,74),-NT(3,25),-NT(7,3));
      CGAL::test_algebraic_structure<NT,Tag>(-NT(5,74),-NT(3,25),-NT(7,3));
      
      CGAL::test_real_embeddable<NT>();
      }
  */
  
  {   // see also  Coercion_traits_test.C
      typedef typename AT::Integer                 I ;
      typedef CGAL::Quotient<typename AT::Integer> QI;
      typedef CGAL::Coercion_traits<I,QI>  CT;
      CGAL_USE_TYPE(CT);
      CGAL_static_assertion((boost::is_same< typename CT::Are_explicit_interoperable,CGAL::Tag_true>::value));
      CGAL_static_assertion((boost::is_same< typename CT::Are_implicit_interoperable,CGAL::Tag_true>::value));
      CGAL_static_assertion((boost::is_same< typename CT::Type,QI>::value));
  }
}

int main() {
#ifdef CGAL_HAS_LEDA_ARITHMETIC_KERNEL
  test_quotient<CGAL::LEDA_arithmetic_kernel>();
#endif
#ifdef CGAL_HAS_CORE_ARITHMETIC_KERNEL
  test_quotient<CGAL::CORE_arithmetic_kernel>();
#endif
#ifdef CGAL_HAS_GMP_ARITHMETIC_KERNEL
  test_quotient<CGAL::GMP_arithmetic_kernel>();
#endif
  return 0;
}
