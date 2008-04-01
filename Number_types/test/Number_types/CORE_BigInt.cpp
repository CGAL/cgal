#include <iostream>

#include <CGAL/basic.h>
#include <cassert>
#ifdef CGAL_USE_CORE

#include <CGAL/CORE_BigInt.h>
#include <CGAL/Test/_test_algebraic_structure.h>
#include <CGAL/Test/_test_real_embeddable.h>
#include <CGAL/Needs_parens_as_product.h>

void test_io(){
    typedef CORE::BigInt NT;
    // MODE ASCII
    {
        std::stringstream ss;
        CGAL::set_ascii_mode(ss);
        ss << CGAL::oformat(NT(1));
        assert( ss.str() == "1");
    }{
        std::stringstream ss;
        CGAL::set_ascii_mode(ss);
        ss << CGAL::oformat(NT(0));
        assert( ss.str() == "0");
    }{
        std::stringstream ss;
        CGAL::set_ascii_mode(ss);
        ss << CGAL::oformat(NT(-1));
        assert( ss.str() == "-1");
    }
    //MODE PRETTY
    {
        std::stringstream ss;
        CGAL::set_pretty_mode(ss);
        ss << CGAL::oformat(NT(1), CGAL::Parens_as_product_tag());
        assert( ss.str() == "1");
    }{
        std::stringstream ss;
        CGAL::set_pretty_mode(ss);
        ss << CGAL::oformat(NT(0),CGAL::Parens_as_product_tag());
        assert( ss.str() == "0");
    }{
        std::stringstream ss;
        CGAL::set_pretty_mode(ss);
        ss << CGAL::oformat(NT(-1), CGAL::Parens_as_product_tag());
        assert( ss.str() == "(-1)");
    }
}

int main() {
    typedef CORE::BigInt NT;
    typedef CGAL::Euclidean_ring_tag Tag;
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
    test_io();
  return 0;
}

#else
int main() { return 0; }
#endif // CGAL_USE_CORE

//EOF
 
