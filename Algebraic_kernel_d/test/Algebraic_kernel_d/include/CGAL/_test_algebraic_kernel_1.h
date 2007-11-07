// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     : Sebastian Limbach <slimbach@mpi-inf.mpg.de>
//                 Michael Hemmer    <hemmer@mpi-inf.mpg.de> 
//
// ============================================================================

#include <CGAL/basic.h>
#include <CGAL/Testsuite/assert.h>
#include <CGAL/Algebraic_kernel_1.h>

#include <CGAL/_test_basic.h>

// Test for the Algebraic_kernel syntax
#ifndef CGAL_TEST_ALGEBRAIC_KERNEL_1_H
#define CGAL_TEST_ALGEBRAIC_KERNEL_1_H

CGAL_BEGIN_NAMESPACE

namespace CGALi {

template< class AK_, class AlgebraicReal1, class Isolator_, class Coefficient_, class Polynomial1, class Boundary_  >
void test_algebraic_kernel_1() {
    typedef AK_            AK;
    typedef AlgebraicReal1 Algebraic_real_1;
    typedef Isolator_      Isolator;
    typedef Coefficient_   Coefficient;
    typedef Polynomial1    Polynomial_1;
    typedef Boundary_      Boundary;
        
    BOOST_STATIC_ASSERT( (::boost::is_same< 
            Algebraic_real_1, typename AK::Algebraic_real_1 >::value) );

    BOOST_STATIC_ASSERT((::boost::is_same<
            Isolator,
            typename AK::Isolator >::value) );
            
    BOOST_STATIC_ASSERT((::boost::is_same< 
            Coefficient, 
            typename AK::Coefficient >::value));
            
    BOOST_STATIC_ASSERT((::boost::is_same<
            Polynomial_1,
            typename AK::Polynomial_1 >::value));
    
    // Test of functors
    // Test AK::Solve_1...
    typename AK::Solve_1 solve_1;
    Polynomial_1 poly1( -4,0,1 );
    Polynomial_1 poly2( 0, 0, 1 );
    std::vector< Algebraic_real_1 > roots_vec;
    std::vector< int > mults_vec;
    
    solve_1( poly1, std::back_inserter( roots_vec ) );        
    CGAL_test_assert( roots_vec.size() == 2 );
    CGAL_test_assert( CGAL::abs( roots_vec[0] ) == CGAL::abs( Algebraic_real_1(2) ) );
    roots_vec.clear();
    
    solve_1( poly1, std::back_inserter( roots_vec ), std::back_inserter( mults_vec ) );
    CGAL_test_assert( roots_vec.size() == 2 );
    CGAL_test_assert( mults_vec.size() == 2 );
    CGAL_test_assert( CGAL::abs( roots_vec[0] ) == CGAL::abs( Algebraic_real_1(2) ) );
    CGAL_test_assert( CGAL::abs( roots_vec[1] ) == CGAL::abs( Algebraic_real_1(2) ) );
    CGAL_test_assert( mults_vec[0] == 1 );
    CGAL_test_assert( mults_vec[1] == 1 );
    roots_vec.clear();
    mults_vec.clear();

    solve_1( poly2, std::back_inserter( roots_vec ), std::back_inserter( mults_vec ) );
    CGAL_test_assert( roots_vec.size() == 1 );
    CGAL_test_assert( mults_vec.size() == 1 );
    CGAL_test_assert( CGAL::abs( roots_vec[0] ) == CGAL::abs( Algebraic_real_1(0) ) );
    CGAL_test_assert( mults_vec[0] == 2 );        
    roots_vec.clear();
    mults_vec.clear();
    
    // Test AK::Sign_at_1
    typename AK::Sign_at_1 sign_at_1;
    typename AK::Polynomial_1 poly4( -2,0,1 );
    solve_1( poly4, std::back_inserter( roots_vec ) );
    typename AK::Polynomial_1 poly3( 0,0,0,1 ); 
    CGAL_test_assert( sign_at_1( poly3, roots_vec[0] ) == CGAL::sign( roots_vec[0] ) );
    CGAL_test_assert( sign_at_1( poly3, roots_vec[1] ) == CGAL::sign( roots_vec[1] ) );
    CGAL_test_assert( sign_at_1( poly3, Algebraic_real_1(0) ) == CGAL::ZERO );  
    roots_vec.clear();
    
    solve_1( poly1, std::back_inserter( roots_vec ) );
    CGAL_test_assert( sign_at_1( poly3, roots_vec[0] ) == CGAL::sign( roots_vec[0] ) );
    CGAL_test_assert( sign_at_1( poly3, roots_vec[1] ) == CGAL::sign( roots_vec[1] ) );
    CGAL_test_assert( sign_at_1( poly3, Algebraic_real_1(0) ) == CGAL::ZERO );  
    roots_vec.clear();
    
    typename AK::Polynomial_1 poly5( 0,0,-1,0,1 );
    typename AK::Algebraic_real_1 algreal1( poly1, Boundary(-3), Boundary(1) );
    typename AK::Algebraic_real_1 algreal2( poly1, Boundary(-1), Boundary(3) );
    CGAL_test_assert( sign_at_1( poly5, algreal2 ) == CGAL::POSITIVE );
    
    
    // Just syntax tests... (TODO)
    // Test AK::Is_square_free_1...
    typename AK::Is_square_free_1 is_square_free_1;
    is_square_free_1( poly1 );
    
    // Test AK::Is_coprime_1...
    typename AK::Is_coprime_1 is_coprime_1;
    is_coprime_1( poly1, poly2 );
        
    // Test AK::Make_square_free_1...
    typename AK::Make_square_free_1 make_square_free_1;
    make_square_free_1( poly1 );
    
    // Test AK::Make_coprime_1...
    typename AK::Make_coprime_1 make_coprime_1;
    Polynomial_1 g, q1, q2;
    make_coprime_1( poly1, poly2, g, q1, q2 );
    
    // Test AK::Square_free_factorization_1...
    typename AK::Square_free_factorization_1 square_free_factorization_1;
    std::vector<Polynomial_1> factors;
    std::vector<int> mults;
    square_free_factorization_1( poly1, std::back_inserter(factors), 
                                     std::back_inserter(mults) );
        
    ////////////////////////////////////////////////////////////////////////////
    
    // (Not only) syntax tests for Algebraic_real_traits
    typedef typename AK::Algebraic_real_traits ART;
    
    BOOST_STATIC_ASSERT((::boost::is_same<
            typename AK::Algebraic_real_1,
            typename ART::Type >::value));
            
    BOOST_STATIC_ASSERT((::boost::is_same<
            Boundary,
            typename ART::Boundary >::value));        
        
    // Create test polynomial
    Polynomial_1 p2( -2,0,1 );
    std::vector< Algebraic_real_1 > roots_vec2;
    
    solve_1( p2, std::back_inserter( roots_vec2 ) );
    
    // Test ART::Boundary_between...
    typename ART::Boundary_between boundary_between;
    CGAL_test_assert( typename ART::Boundary( -2 ) < boundary_between( roots_vec2[0], roots_vec2[1] ) );
    CGAL_test_assert( typename ART::Boundary(  2 ) > boundary_between( roots_vec2[0], roots_vec2[1] ) );
    
    // Test ART::Lower_boundary
    typename ART::Lower_boundary lower_boundary;
    CGAL_test_assert( lower_boundary( roots_vec2[0] ) < typename ART::Boundary(-1) );
    CGAL_test_assert( lower_boundary( roots_vec2[1] ) < typename ART::Boundary( 2) );

    // Test ART::Upper_boundary
    typename ART::Upper_boundary upper_boundary;
    CGAL_test_assert( upper_boundary( roots_vec2[0] ) > typename ART::Boundary(-1) );
    CGAL_test_assert( upper_boundary( roots_vec2[1] ) > typename ART::Boundary( 1) );
    
    // Test ART::Refine
    typename ART::Refine refine;
    Algebraic_real_1 ar = roots_vec2[1];
    typename ART::Boundary old_lower_boundary = ar.low();
    typename ART::Boundary old_upper_boundary = ar.high(); 

    refine( ar );
    
    CGAL_test_assert( old_lower_boundary <= lower_boundary( ar ) );
    CGAL_test_assert( old_upper_boundary >= upper_boundary( ar ) );
    typename ART::Boundary interval_size_old = CGAL::abs( old_upper_boundary - old_lower_boundary );
    typename ART::Boundary interval_size_new = CGAL::abs( upper_boundary( ar ) - lower_boundary( ar ) );
    CGAL_test_assert( interval_size_new * typename ART::Boundary(2) <= interval_size_old );

    refine( ar, 100 );
    CGAL_test_assert( CGAL::abs( upper_boundary( ar ) - lower_boundary( ar ) ) < 
        (typename ART::Boundary(1) / POLYNOMIAL::ipower(typename ART::Boundary(2), 99 )) );

}

} //namespace CGALi

CGAL_END_NAMESPACE

#endif //CGAL_TEST_ALGEBRAIC_KERNEL_1_H
