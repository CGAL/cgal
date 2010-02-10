

#include <CGAL/basic.h>

#ifdef CGAL_USE_GMP

#include <iostream>
#include <CGAL/Gmpz.h>
#include <CGAL/Test/_test_algebraic_structure.h>
#include <CGAL/Test/_test_real_embeddable.h>

int main() {
{   typedef CGAL::Gmpz NT;
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

    //test additional functions
    NT a(4);
    assert((a >>= 1) == 2); 
    assert((a <<= 1) == 4); 

    assert( a++ == 4); 
    assert( a   == 5); 
    assert( ++a == 6);

    assert( a-- == 6); 
    assert( a   == 5); 
    assert( --a == 4);

    assert( (NT(5) & NT(4)) == NT(4));
    assert( (NT(5) | NT(4)) == NT(5));
    assert( (NT(5) ^ NT(4)) == NT(1));
}


  return 0;
}

#else 
int main()
{
  return 0;
}
#endif //CGAL_USE_GMP
