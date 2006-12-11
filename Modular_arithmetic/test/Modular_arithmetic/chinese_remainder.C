// Author(s)     : Michael Hemmer <mhemmer@uni-mainz.de>


#include <CGAL/basic.h>

#include <CGAL/Testsuite/assert.h>
#include <CGAL/chinese_remainder.h>
#include <cstdlib>

#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
#endif

#ifdef CGAL_USE_CORE
#include <CGAL/CORE_BigInt.h>
#endif

template <class NT>
void test_chinese_remainder(NT x){
    typedef CGAL::Algebraic_structure_traits<NT> AST;
    typename AST::Mod mod;
    
    NT m1 = 23;
    NT m2 = 17;
    NT m3 = 29;
    NT u1 = mod(x,m1);
    NT u2 = mod(x,m2);
    NT u3 = mod(x,m3);
    NT m,u;

    CGAL::chinese_remainder(m1,u1,m2,u2,m,u);
    CGAL::chinese_remainder(m ,u ,m3,u3,m,u);

    CGAL_test_assert( m  == m1*m2*m3 );
    CGAL_test_assert( x  == u );
}

template <class NT>
void test_chinese_remainder(){
    test_chinese_remainder(NT(0));
    test_chinese_remainder(NT(1));
    test_chinese_remainder(NT(-1));
    test_chinese_remainder(NT(23));
    test_chinese_remainder(NT(17));
    test_chinese_remainder(NT(-29));
    test_chinese_remainder(NT(2456));
    test_chinese_remainder(NT(-2456));
}

int main(){
    test_chinese_remainder<int>();

#ifdef CGAL_USE_LEDA
    test_chinese_remainder<leda_integer>();
#endif

#ifdef CGAL_USE_CORE
    test_chinese_remainder<CORE::BigInt>();
#endif

    return 0; 
}
