#include <iostream>

#include <CGAL/config.h>
#ifdef CGAL_USE_CORE

#include <CGAL/CORE_BigFloat.h>
#include <CGAL/Test/_test_algebraic_structure.h>
#include <CGAL/Test/_test_real_embeddable.h>

int main() {
    typedef CORE::BigFloat NT;
    typedef CGAL::Field_with_kth_root_tag Tag;
    typedef CGAL::Tag_false Is_exact;
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

    {
      CORE::BigInt m("-3295371043742771943824");
      CORE::BigFloat bf(m,16,-5);
      CGAL::to_interval(bf);
    }

  return 0;
}

#else // CGAL_USE_CORE
int main() { return 0; }
#endif // CGAL_USE_CORE

//EOF

