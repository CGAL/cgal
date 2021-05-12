
#include <CGAL/config.h>

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
    assert((a >> 1) == 2);
    assert( a == 4);
    assert((a << 1) == 8);
    assert( a == 4);
    assert(((a >>= 1) >>= 1)== 1);
    assert( a == 1);
    assert(((a <<= 1) <<= 1) == 4);
    assert( a == 4);

    assert( a++ == 4);
    assert( a   == 5);
    assert( ++(++a) == 7);
    assert( a == 7);

    assert( a-- == 7);
    assert( a   == 6);
    assert( --(--a) == 4);

    NT b(5);
    assert( (b & a) == 4);
    assert( (b | a) == 5);
    assert( (b ^ a) == 1);
    assert( a==4 );
    assert( b==5 );

    assert( (b &= NT(4)) == 4);
    assert( b == 4);
    assert( (b |= NT(5)) == 5);
    assert( b == 5);
    assert( (b ^= NT(4)) == 1);
    assert( b==1 );
}


  return 0;
}

#else
int main()
{
  return 0;
}
#endif //CGAL_USE_GMP
