#include <iostream>

#include <CGAL/basic.h>
#ifdef CGAL_USE_LEDA
#include <CGAL/leda_rational.h>
#include <CGAL/Nef_polynomial.h>

#include <CGAL/_test_algebraic_structure.h>
#include <CGAL/_test_real_embeddable.h>

int main() {
    typedef leda_rational CT;
    typedef CGAL::Nef_polynomial<CT> NT;
    typedef CGAL::Integral_domain_without_division_tag Tag;
    //TODO    typedef CGAL::Unique_factorization_domain_tag Tag;
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
  
    CGAL::test_real_embeddable<NT>();

//     // backward compatiblity
//     typedef CGAL::Number_type_traits<NT> NTT;
//     CGAL_test_assert(  CGAL::check_tag(NTT::Has_division()) );
//     CGAL_test_assert(  CGAL::check_tag(NTT::Has_gcd()) );
//     CGAL_test_assert(  CGAL::check_tag(NTT::Has_sqrt()) );
//     CGAL_test_assert(  CGAL::check_tag(NTT::Has_exact_ring_operations()) );
//     CGAL_test_assert(! CGAL::check_tag(NTT::Has_exact_division()) );
//     CGAL_test_assert(! CGAL::check_tag(NTT::Has_exact_sqrt()) );
    
  return 0;
}

#else
int main() { return 0; }
#endif
