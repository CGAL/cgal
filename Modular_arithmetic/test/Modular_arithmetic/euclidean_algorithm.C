// Author(s)     : Michael Hemmer <mhemmer@uni-mainz.de>


#include <CGAL/basic.h>

#include <CGAL/Testsuite/assert.h>
#include <CGAL/euclidean_algorithm.h>

#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
#endif

#ifdef CGAL_USE_CORE
#include <CGAL/CORE_BigInt.h>
#endif

#include <cstdlib>


template <class NT>
void test_euclidean_algorithm(){
    CGAL_test_assert( CGAL::euclidean_algorithm(NT( 0),NT(0)) == NT( 0));
    CGAL_test_assert( CGAL::euclidean_algorithm(NT( 7),NT(0)) == NT( 7));
    CGAL_test_assert( CGAL::euclidean_algorithm(NT(-7),NT(0)) == NT( 7));
    CGAL_test_assert( CGAL::euclidean_algorithm(NT( 0),NT(7)) == NT( 7));
    CGAL_test_assert( CGAL::euclidean_algorithm(NT(0),NT(-7)) == NT( 7));

    CGAL_test_assert( CGAL::euclidean_algorithm(NT( 1),NT(1)) == NT( 1));
    CGAL_test_assert( CGAL::euclidean_algorithm(NT( 7),NT(1)) == NT( 1));
    CGAL_test_assert( CGAL::euclidean_algorithm(NT(-7),NT(1)) == NT( 1));
    CGAL_test_assert( CGAL::euclidean_algorithm(NT( 1),NT(7)) == NT( 1));
    CGAL_test_assert( CGAL::euclidean_algorithm(NT(1),NT(-7)) == NT( 1));

    CGAL_test_assert( CGAL::euclidean_algorithm(NT( 7),NT(7)) == NT( 7));

    CGAL_test_assert( CGAL::euclidean_algorithm(NT( 3),NT(1)) == NT( 1));
    CGAL_test_assert( CGAL::euclidean_algorithm(NT( 3),NT(6)) == NT( 3));
    CGAL_test_assert( CGAL::euclidean_algorithm(NT( 6),NT(9)) == NT( 3));
    CGAL_test_assert( CGAL::euclidean_algorithm(NT( 9),NT(15)) == NT( 3));
    CGAL_test_assert( CGAL::euclidean_algorithm(NT( 15),NT(24)) == NT( 3));
    CGAL_test_assert( CGAL::euclidean_algorithm(NT( 24),NT(39)) == NT( 3));
    CGAL_test_assert( CGAL::euclidean_algorithm(NT( 39),NT(63)) == NT( 3));
    CGAL_test_assert( CGAL::euclidean_algorithm(NT( 6),NT(3)) == NT( 3));
    CGAL_test_assert( CGAL::euclidean_algorithm(NT( 9),NT(6)) == NT( 3));
    CGAL_test_assert( CGAL::euclidean_algorithm(NT( 15),NT(9)) == NT( 3));
    CGAL_test_assert( CGAL::euclidean_algorithm(NT( 24),NT(15)) == NT( 3));
    CGAL_test_assert( CGAL::euclidean_algorithm(NT( 39),NT(24)) == NT( 3));
    CGAL_test_assert( CGAL::euclidean_algorithm(NT( 63),NT(39)) == NT( 3));

    CGAL_test_assert( CGAL::euclidean_algorithm(NT( -3),NT(6)) == NT( 3));
    CGAL_test_assert( CGAL::euclidean_algorithm(NT( -6),NT(9)) == NT( 3));
    CGAL_test_assert( CGAL::euclidean_algorithm(NT( -9),NT(15)) == NT( 3));
    CGAL_test_assert( CGAL::euclidean_algorithm(NT( -15),NT(24)) == NT( 3));
    CGAL_test_assert( CGAL::euclidean_algorithm(NT( -6),NT(3)) == NT( 3));
    CGAL_test_assert( CGAL::euclidean_algorithm(NT( -9),NT(6)) == NT( 3));
    CGAL_test_assert( CGAL::euclidean_algorithm(NT( -15),NT(9)) == NT( 3));
    CGAL_test_assert( CGAL::euclidean_algorithm(NT( -24),NT(15)) == NT( 3));

    CGAL_test_assert( CGAL::euclidean_algorithm(NT( 3),NT(-1)) == NT( 1));
    CGAL_test_assert( CGAL::euclidean_algorithm(NT( 3),NT(-6)) == NT( 3));
    CGAL_test_assert( CGAL::euclidean_algorithm(NT( 6),NT(-9)) == NT( 3));
    CGAL_test_assert( CGAL::euclidean_algorithm(NT( 6),NT(-3)) == NT( 3));
    CGAL_test_assert( CGAL::euclidean_algorithm(NT( 9),NT(-6)) == NT( 3));

    CGAL_test_assert( CGAL::euclidean_algorithm(NT( -3),NT(-6)) == NT( 3));
    CGAL_test_assert( CGAL::euclidean_algorithm(NT( -6),NT(-9)) == NT( 3));
    CGAL_test_assert( CGAL::euclidean_algorithm(NT( -6),NT(-3)) == NT( 3));
    CGAL_test_assert( CGAL::euclidean_algorithm(NT( -9),NT(-6)) == NT( 3));
}


int main(){
    test_euclidean_algorithm<int>();
    test_euclidean_algorithm<leda_integer>();
    test_euclidean_algorithm<CORE::BigInt>();
    return 0; 
}
