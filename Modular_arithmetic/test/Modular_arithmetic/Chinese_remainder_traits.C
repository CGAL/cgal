// Author(s)     : Michael Hemmer <mhemmer@uni-mainz.de>

#include <CGAL/basic.h>
#include <CGAL/Testsuite/assert.h>
#include <CGAL/Chinese_remainder_traits.h>

#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
#endif

#ifdef CGAL_USE_CORE
#include <CGAL/CORE_BigInt.h>
#endif

template <class CRT>
void test_chinese_remainder_traits(){
    typedef typename CRT::T T;
    typedef typename CRT::Scalar_type Scalar_type;
    typedef typename CRT::Chinese_remainder Chinese_remainder;
    Chinese_remainder chinese_remainder;
    
    T x(-574);
    Scalar_type m1 = Scalar_type(23);
    Scalar_type m2 = Scalar_type(17);
    Scalar_type m3 = Scalar_type(29);
    T u1 = -T(22); 
    T u2 = -T(13); 
    T u3 = -T(23); 
    Scalar_type m;
    T u;

    chinese_remainder(m1,u1,m2,u2,m,u);
    chinese_remainder(m ,u ,m3,u3,m,u);

    CGAL_test_assert( m  == m1*m2*m3 );
    CGAL_test_assert( x  == u );

}


int main(){
  
    test_chinese_remainder_traits<CGAL::Chinese_remainder_traits<int> >();
    typedef CGAL::Sqrt_extension<int,int>       Extn_1;
    typedef CGAL::Sqrt_extension<Extn_1,int>    Extn_2;
    typedef CGAL::Sqrt_extension<Extn_1,Extn_1> Extn_n2;
        
    test_chinese_remainder_traits<CGAL::Chinese_remainder_traits<Extn_1 > >();
    test_chinese_remainder_traits<CGAL::Chinese_remainder_traits<Extn_2 > >();
    test_chinese_remainder_traits<CGAL::Chinese_remainder_traits<Extn_n2 > >();

    typedef CGAL::Polynomial<int> Poly_1;
    typedef CGAL::Polynomial<Poly_1> Poly_2;
    test_chinese_remainder_traits<CGAL::Chinese_remainder_traits<Poly_1 > >();
    test_chinese_remainder_traits<CGAL::Chinese_remainder_traits<Poly_2 > >();
    

#ifdef CGAL_USE_CORE
    test_chinese_remainder_traits<
    CGAL::Chinese_remainder_traits<CORE::BigInt> >();
#endif

#ifdef CGAL_USE_LEDA
    test_chinese_remainder_traits<
    CGAL::Chinese_remainder_traits<leda::integer> >();
#endif


}
