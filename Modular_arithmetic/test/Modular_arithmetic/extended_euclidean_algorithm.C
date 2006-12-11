// Author(s)     : Michael Hemmer <mhemmer@uni-mainz.de>


#include <CGAL/basic.h>

#include <CGAL/Testsuite/assert.h>
#include <CGAL/extended_euclidean_algorithm.h>
#include <cstdlib>


template <class NT>
void test_extended_euclidean_algorithm
(const NT& a, const NT& b, const NT& g_){
    NT u,v; 
    NT g = CGAL::extended_euclidean_algorithm(a,b,u,v);
    CGAL_test_assert(g_ == g) ;
    CGAL_test_assert(g_ == u*a+v*b);
}


template <class NT>
void test_extended_euclidean_algorithm(){ 

    test_extended_euclidean_algorithm(NT( 0),NT(0) , NT( 0));
    test_extended_euclidean_algorithm(NT( 7),NT(0) , NT( 7));
    test_extended_euclidean_algorithm(NT(-7),NT(0) , NT( 7));
    test_extended_euclidean_algorithm(NT( 0),NT(7) , NT( 7));
    test_extended_euclidean_algorithm(NT(0),NT(-7) , NT( 7));

    test_extended_euclidean_algorithm(NT( 1),NT(1) , NT( 1));
    test_extended_euclidean_algorithm(NT( 7),NT(1) , NT( 1));
    test_extended_euclidean_algorithm(NT(-7),NT(1) , NT( 1));
    test_extended_euclidean_algorithm(NT( 1),NT(7) , NT( 1));
    test_extended_euclidean_algorithm(NT(1),NT(-7) , NT( 1));
    
    test_extended_euclidean_algorithm(NT( 7),NT(7) , NT( 7));
    
    test_extended_euclidean_algorithm(NT( 3),NT(1) , NT( 1));
    test_extended_euclidean_algorithm(NT( 3),NT(6) , NT( 3));
    test_extended_euclidean_algorithm(NT( 6),NT(9) , NT( 3));
    test_extended_euclidean_algorithm(NT( 9),NT(15) , NT( 3));
    test_extended_euclidean_algorithm(NT( 15),NT(24) , NT( 3));
    test_extended_euclidean_algorithm(NT( 24),NT(39) , NT( 3));
    test_extended_euclidean_algorithm(NT( 39),NT(63) , NT( 3));
    test_extended_euclidean_algorithm(NT( 6),NT(3) , NT( 3));
    test_extended_euclidean_algorithm(NT( 9),NT(6) , NT( 3));
    test_extended_euclidean_algorithm(NT( 15),NT(9) , NT( 3));
    test_extended_euclidean_algorithm(NT( 24),NT(15) , NT( 3));
    test_extended_euclidean_algorithm(NT( 39),NT(24) , NT( 3));
    test_extended_euclidean_algorithm(NT( 63),NT(39) , NT( 3));
    
    test_extended_euclidean_algorithm(NT( -3),NT(6) , NT( 3));
    test_extended_euclidean_algorithm(NT( -6),NT(9) , NT( 3));
    test_extended_euclidean_algorithm(NT( -9),NT(15) , NT( 3));
    test_extended_euclidean_algorithm(NT( -15),NT(24) , NT( 3));
    test_extended_euclidean_algorithm(NT( -6),NT(3) , NT( 3));
    test_extended_euclidean_algorithm(NT( -9),NT(6) , NT( 3));
    test_extended_euclidean_algorithm(NT( -15),NT(9) , NT( 3));
    test_extended_euclidean_algorithm(NT( -24),NT(15) , NT( 3));

    test_extended_euclidean_algorithm(NT( 3),NT(-1) , NT( 1));
    test_extended_euclidean_algorithm(NT( 3),NT(-6) , NT( 3));
    test_extended_euclidean_algorithm(NT( 6),NT(-9) , NT( 3));
    test_extended_euclidean_algorithm(NT( 6),NT(-3) , NT( 3));
    test_extended_euclidean_algorithm(NT( 9),NT(-6) , NT( 3));

    test_extended_euclidean_algorithm(NT( -3),NT(-6) , NT( 3));
    test_extended_euclidean_algorithm(NT( -6),NT(-9) , NT( 3));
    test_extended_euclidean_algorithm(NT( -6),NT(-3) , NT( 3));
    test_extended_euclidean_algorithm(NT( -9),NT(-6), NT( 3));
}


int main(){
    test_extended_euclidean_algorithm<int>();
    return 0; 
}
