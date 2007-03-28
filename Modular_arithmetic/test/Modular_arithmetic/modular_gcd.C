// Author(s)     : Michael Hemmer <mhemmer@uni-mainz.de>

/*! \file CGAL/Modular.C
  test for number type modul 
*/

#include <CGAL/basic.h>
#include <CGAL/Testsuite/assert.h>
#include <CGAL/Modular.h>
#include <CGAL/modular_gcd.h>
#include <CGAL/Polynomial.h>

#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
#endif // CGAL_USE_LEDA

#include <cstdlib>

#include <boost/type_traits.hpp>

int main()
{  
    { 
        typedef leda::integer Integer;
        typedef CGAL::Polynomial<Integer> Polynomial;
        Polynomial p1(123,431,2134);
        Polynomial p2(123,421,234);
        Polynomial g(1);
        Polynomial result = modular_gcd_utcf(p1,p2);
        //std::cout <<" result    : " << result <<std::endl;
        //std::cout <<" true gcd  : " << g <<std::endl;
        CGAL_test_assert(result == g);
        
    }
    {  // unlucky prime test 
        typedef leda::integer Integer;
        typedef CGAL::Polynomial<Integer> Polynomial;
        
        Polynomial f1(5,234,445);
        Polynomial f2(12,-234,345);
        f1 *= Polynomial(Integer(CGAL::primes[0]+3),Integer(1));
        f2 *= Polynomial(Integer(3),Integer(1));
        Polynomial g(13,96,2345);
        Polynomial p1 = f1*g;
        Polynomial p2 = f2*g;
        
        Polynomial result = modular_gcd_utcf(p1,p2);
        //std::cout <<" result    : " << result <<std::endl;
        //std::cout <<" true gcd  : " << g <<std::endl;
        CGAL_test_assert(result == g);
        
    }
 
    {
        typedef leda::integer Integer;
        typedef CGAL::Polynomial<Integer> Polynomial;
        
        Polynomial f1(5,234,-26,243,745);
        Polynomial f2(12,-234,26,243,-731);
        Polynomial g(13,-5676,234,96);
        Polynomial p1 = Polynomial(8)*f1*f1*g;
        Polynomial p2 = Polynomial(5)*f2*f2*g;
        
        Polynomial result = modular_gcd_utcf(p1,p2);
        //std::cout <<" result    : " << result <<std::endl;
        //std::cout <<" true gcd  : " << g <<std::endl;
        CGAL_test_assert(result == g);
        
    }
    {   // test for sqrt 
        typedef leda::integer Integer;
        typedef CGAL::Sqrt_extension<Integer,Integer> EXT;
        typedef CGAL::Polynomial<EXT> Polynomial;
        
        Integer root(Integer(789234)); 
        Polynomial f1(EXT(235143,-2234,root),EXT(232543,-2334,root),EXT(235403,-2394,root),EXT(235483,-2364,root),EXT(223443,-2234,root));
        Polynomial f2(EXT(25143,-2134,root),EXT(212543,-2315,root),EXT(255453,-5394,root),EXT(535483,-2354,root),EXT(22333,-2214,root));
        Polynomial g(EXT(215143,-2134,root),EXT(2122422543,-2115,root),EXT(255453,-1394,root),EXT(135483,-2354,root),EXT(7));
        g=g*g;g=g*g;g=g*g;
        
        Polynomial p1 = Polynomial(8)*f1*f1*g;
        Polynomial p2 = Polynomial(5)*f2*f2*g;
        
        Polynomial result = modular_gcd_utcf(p1,p2);
        //std::cout <<" result    : " << result <<std::endl;
        //std::cout <<" true gcd  : " << g <<std::endl;
        CGAL_test_assert(result == g);
    }

    {   // test for sqrt / denom for algebraic integer
        CGAL::set_pretty_mode(std::cout);
        typedef leda::integer Integer;
        typedef CGAL::Sqrt_extension<Integer,Integer> EXT;
        typedef CGAL::Polynomial<EXT> Polynomial;
        
        Integer root(4*3); 
        Polynomial f1(EXT(0,1,root),EXT(2));
        Polynomial f2(EXT(0,1,root),EXT(4));
        Polynomial g(EXT(0,-1,root),EXT(2));
        

        //std::cout <<" f1        : " << f1 <<std::endl;
        //std::cout <<" f2        : " << f2 <<std::endl;
        //std::cout <<" g         : " << g <<std::endl;

        Polynomial p1 = f1*g;
        Polynomial p2 = f2*g;
        
        
        //std::cout <<" p1        : " << p1 <<std::endl;
        //std::cout <<" p2        : " << p2 <<std::endl;
        //std::cout << std::endl;
        
        Polynomial result = modular_gcd_utcf(p1,p2);
        //std::cout <<" result    : " << result <<std::endl;
        //std::cout <<" true gcd  : " << g <<std::endl;
        
        CGAL_test_assert(result == g);
    } 
    




    return 0 ;
}


