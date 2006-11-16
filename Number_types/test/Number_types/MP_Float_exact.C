#ifdef CGAL_MP_FLOAT_ALLOW_INEXACT
#undef CGAL_MP_FLOAT_ALLOW_INEXACT
#endif 

#include <CGAL/basic.h>


#include <iostream>
#include <CGAL/MP_Float.h>
#include <CGAL/_test_algebraic_structure.h>
#include <CGAL/_test_real_embeddable.h>

int main() {
    typedef CGAL::MP_Float NT;

    typedef CGAL::Euclidean_ring_tag Tag;
    typedef CGAL::Tag_true Is_exact;

    BOOST_STATIC_ASSERT( CGAL::CGALi::Is_integral_domain_without_division<NT>::value);
    BOOST_STATIC_ASSERT( CGAL::CGALi::Is_integral_domain<NT>::value);
    BOOST_STATIC_ASSERT( CGAL::CGALi::Is_unique_factorization_domain<NT>::value);
    BOOST_STATIC_ASSERT( CGAL::CGALi::Is_euclidean_ring<NT>::value);
    BOOST_STATIC_ASSERT(!CGAL::CGALi::Is_field<NT>::value);
    BOOST_STATIC_ASSERT(!CGAL::CGALi::Is_field_with_sqrt<NT>::value);
    BOOST_STATIC_ASSERT(!CGAL::CGALi::Is_field_with_root_of<NT>::value);

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

   
  return 0;
}
