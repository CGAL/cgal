//#define CGAL_MP_FLOAT_ALLOW_INEXACT // TODO: How can I test both cases in one 
                                    //         file?
#include <iostream>
#include <CGAL/MP_Float.h>
#include <CGAL/_test_algebraic_structure.h>
#include <CGAL/_test_real_embeddable.h>

int main() {
    typedef CGAL::MP_Float NT;

#ifdef CGAL_MP_FLOAT_ALLOW_INEXACT
    typedef CGAL::Field_with_sqrt_tag Tag;
    typedef CGAL::Tag_false Is_exact;
#else
    typedef CGAL::Euclidean_ring_tag Tag;
    typedef CGAL::Tag_true Is_exact;
#endif

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
