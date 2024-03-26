#include <iostream>
#include <sstream>

#include <CGAL/config.h>
#include <cassert>
#ifdef CGAL_USE_CORE

#include <CGAL/CORE_BigRat.h>
#include <CGAL/Needs_parens_as_product.h>
#include <CGAL/Test/_test_algebraic_structure.h>
#include <CGAL/Test/_test_real_embeddable.h>
#include <CGAL/Test/_test_fraction_traits.h>
#include <CGAL/Test/_test_rational_traits.h>

void test_io(){
    typedef CORE::BigRat NT;
    // MODE ASCII
    {
        std::stringstream ss;
        CGAL::IO::set_ascii_mode(ss);
        ss << CGAL::IO::oformat(NT(1));
        std::cout << ss.str()<<std::endl;
        assert( ss.str() == "1");
    }{
        std::stringstream ss;
        CGAL::IO::set_ascii_mode(ss);
        ss << CGAL::IO::oformat(NT(0));
        assert( ss.str() == "0");
    }{
        std::stringstream ss;
        CGAL::IO::set_ascii_mode(ss);
        ss << CGAL::IO::oformat(NT(-1));
        assert( ss.str() == "-1");
    }
    //MODE PRETTY
    {
        std::stringstream ss;
        CGAL::IO::set_pretty_mode(ss);
        ss << CGAL::IO::oformat(NT(2), CGAL::Parens_as_product_tag());
        std::cout << "|" << ss.str() << "|" << std::endl;
        assert( ss.str() == "2");
    }{
        std::stringstream ss;
        CGAL::IO::set_pretty_mode(ss);
        ss << CGAL::IO::oformat(NT(1)/NT(2),CGAL::Parens_as_product_tag());
        assert( ss.str() == "(1/2)");
    }{
        std::stringstream ss;
        CGAL::IO::set_pretty_mode(ss);
        ss << CGAL::IO::oformat(NT(-2), CGAL::Parens_as_product_tag());
        assert( ss.str() == "(-2)");
    }
}

int main() {
    typedef CORE::BigRat NT;
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

    CGAL::test_real_embeddable<NT>();

    CGAL::test_fraction_traits<NT>();
    // backward compatibility
    CGAL::test_rational_traits<NT>();

    test_io();

  return 0;
}

#else
int main() { return 0; }
#endif // CGAL_USE_CORE

//EOF

