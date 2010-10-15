
#include <CGAL/basic.h>

#ifdef CGAL_USE_GMP
#include <iostream>
#include <CGAL/Gmpq.h>
#include <CGAL/Test/_test_algebraic_structure.h>
#include <CGAL/Test/_test_real_embeddable.h>
#include <CGAL/Test/_test_fraction_traits.h>
#include <CGAL/Test/_test_rational_traits.h>
#include <CGAL/to_rational.h>

int main() {
    {   
        typedef CGAL::Gmpq NT;
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
        // backward compatiblity
        CGAL::test_rational_traits<NT>();      
    }

    return 0;
}

#else 
int main()
{
  return 0;
}
#endif //CGAL_USE_GMP
