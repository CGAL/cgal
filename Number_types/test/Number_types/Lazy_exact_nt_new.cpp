#include <iostream>

#include <CGAL/basic.h>
#include <cassert>
#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/Test/_test_algebraic_structure.h>
#include <CGAL/Test/_test_real_embeddable.h>
#include <CGAL/Test/_test_fraction_traits.h>
#include <CGAL/Test/_test_rational_traits.h>
#include <CGAL/number_utils.h>
#ifdef CGAL_USE_CORE
#include <CGAL/CORE_Expr.h>
#endif

template< class AK >
void test_lazy_exact_nt() {
    {
        typedef typename AK::Field_with_sqrt ET;
        typedef CGAL::Algebraic_structure_traits< ET > AST;
        typedef typename AST::Algebraic_category Tag;
        typedef CGAL::Lazy_exact_nt< ET > NT;
        typedef typename AST::Is_exact Is_exact;

        CGAL::test_algebraic_structure<NT,Tag, Is_exact>();
        CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(4),NT(6),NT(15));
        CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(-4),NT(6),NT(15));
        CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(4),NT(-6),NT(15));
        CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(-4),NT(-6),NT(15));
        CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(4),NT(6),NT(-15));
        CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(-4),NT(6), NT(15));
        CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(4),NT(-6),NT(-15));
        CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(-4),NT(-6),NT(-15));
        
        CGAL::test_real_embeddable<NT>();
    }{
        typedef typename AK::Rational ET;

        typedef CGAL::Algebraic_structure_traits< ET > AST;
        typedef typename AST::Algebraic_category Tag;
        typedef CGAL::Lazy_exact_nt< ET > NT;
        typedef typename AST::Is_exact Is_exact;

        CGAL::test_algebraic_structure<NT,Tag, Is_exact>();
        CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(4),NT(6),NT(15));
        CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(-4),NT(6),NT(15));
        CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(4),NT(-6),NT(15));
        CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(-4),NT(-6),NT(15));
        CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(4),NT(6),NT(-15));
        CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(-4),NT(6), NT(15));
        CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(4),NT(-6),NT(-15));
        CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(-4),NT(-6),NT(-15));
        
        CGAL::test_real_embeddable<NT>();

        CGAL::test_fraction_traits<NT>();
        CGAL::test_rational_traits<NT>();
    }
    {        
        typedef typename AK::Integer ET;

        typedef CGAL::Algebraic_structure_traits< ET > AST;
        typedef typename AST::Algebraic_category Tag;
        typedef CGAL::Lazy_exact_nt< ET > NT;
        typedef typename AST::Is_exact Is_exact;

        CGAL::test_algebraic_structure<NT,Tag, Is_exact>();
        CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(4),NT(6),NT(15));
        CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(-4),NT(6),NT(15));
        CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(4),NT(-6),NT(15));
        CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(-4),NT(-6),NT(15));
        CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(4),NT(6),NT(-15));
        CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(-4),NT(6), NT(15));
        CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(4),NT(-6),NT(-15));
        CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(-4),NT(-6),NT(-15));
        
        CGAL::test_real_embeddable<NT>();
    }
    
    {   // see also  Coercion_traits_test.C
        typedef CGAL::Lazy_exact_nt< typename AK::Integer > LI;
        typedef CGAL::Lazy_exact_nt< typename AK::Rational > LR;
        typedef CGAL::Coercion_traits<LI,LR> CT;
        CGAL_static_assertion((boost::is_same< typename CT::Are_implicit_interoperable,CGAL::Tag_true>::value));
        CGAL_static_assertion((boost::is_same< typename CT::Are_explicit_interoperable,CGAL::Tag_true>::value));
        CGAL_static_assertion((boost::is_same< typename CT::Type,LR>::value));
        
        LI  i(4);
        LR  r(4);
        typename CT::Cast cast;
        assert( cast ( (i*i+i) / i-i ) == LR(1));
        assert( cast ( (i*i+r) / i-i ) == LR(1));
        assert( cast ( (i*r+r) / i-i ) == LI(1));
    }{  // see also  Coercion_traits_test.C
#ifdef CGAL_USE_LEDA
#ifdef CGAL_USE_CORE
        typedef CGAL::Lazy_exact_nt<leda_integer  > T1;
        typedef CGAL::Lazy_exact_nt<CORE::Expr    > T2;
        typedef CGAL::Coercion_traits<T1,T2> CT;
        CGAL_static_assertion((boost::is_same< typename CT::Are_implicit_interoperable,CGAL::Tag_false>::value));
        CGAL_static_assertion((boost::is_same< typename CT::Are_explicit_interoperable,CGAL::Tag_false>::value));
#endif
#endif
    }  
}

int main() {
#ifdef CGAL_USE_LEDA
  test_lazy_exact_nt< CGAL::LEDA_arithmetic_kernel >();
#endif
#ifdef CGAL_USE_CORE
  test_lazy_exact_nt< CGAL::CORE_arithmetic_kernel >();
#endif
  return 0;
}

