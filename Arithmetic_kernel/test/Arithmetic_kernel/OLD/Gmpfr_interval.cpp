#include <iostream>
#include <CGAL/basic.h>
#ifdef CGAL_HAS_GMP_ARITHMETIC_KERNEL
#include <CGAL/Gmpfr_interval.h>
#include <CGAL/Test/_test_algebraic_structure.h>
#include <CGAL/Test/_test_real_embeddable.h>
#include <CGAL/Test/_test_interval.h>

int main() {
  std::cout << "TEST Gmpfr_interval "; std::cout.flush();
    typedef CGAL::Gmpfr_interval NT;
    typedef CGAL::Field_with_sqrt_tag Tag;
    typedef CGAL::Tag_false           Is_exact;
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
    CGAL::test_interval<NT>();
    
    std::cout << "OK" << std::endl;
    return 0;
}

#else
int main() { return 0; }
#endif
