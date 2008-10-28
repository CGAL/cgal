// ============================================================================
//
// Copyright (c) 2001-2006 Max-Planck-Institut Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of EXACUS (http://www.mpi-inf.mpg.de/projects/EXACUS/);
// you may redistribute it under the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with EXACUS.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// ----------------------------------------------------------------------------
//
// Library       : CGAL
// File          : test/modular_gcd.C
// CGAL_release   : $Name:  $
// Revision      : $Revision$
// Revision_date : $Date$
//
// Author(s)     : Dominik Huelse <dominik.huelse@gmx.de>
//                 Michael Hemmer <hemmer@mpi-inf.mpg.de>
//
// ============================================================================


#include <CGAL/basic.h>

#include <vector>
#include <cassert>

#include <CGAL/Polynomial/modular_gcd.h>
#include <CGAL/gen_polynomials.h>
#include <CGAL/Random.h>

#include <CGAL/Arithmetic_kernel.h>   


//#define WITH_OUTPUT 1  

// test algorithm modular_gcd
CGAL::Random my_random(4711);

template<class NT>
void gcd_utcf_test(const NT& f, const NT& g, const NT& d) {
    NT tmp = CGAL::CGALi::modular_gcd_utcf_algorithm_M(f, g);
    //     NT tmp = CGAL::gcd_utcf(f, g);
#ifdef WITH_OUTPUT
    std::cout << "\nf(x) = " << f;
    std::cout << "\ng(x) = " << g;
    std::cout << "\ngcd_utcf(f,g) = " << tmp;
    std::cout << "\nd        = " << CGAL::canonicalize(d) << "\n";
#endif
    if(tmp != NT(-1))    // NT(-1) when primes exhausted
        assert( 
                CGAL::canonicalize(tmp) 
                == 
                CGAL::canonicalize(d) );
}

template<class NT>
void sqrt_gcd_utcf_test(const NT& f, const NT& g, const NT& d) {
//    NT tmp = CGAL::CGALi::modular_gcd_utcf_dfai(f, g);
    NT tmp = CGAL::CGALi::modular_gcd_utcf_with_wang(f, g);
 
#ifdef WITH_OUTPUT
    std::cout << "\nf(x) = " << f;  
    std::cout << "\ng(x) = " << g;
    std::cout << "\ngcd_utcf(f,g) = " << tmp;
    std::cout << "\nd        = " << CGAL::canonicalize(d) << "\n";
#endif
    if(tmp != NT(-1))    // NT(-1) when primes exhausted
        assert( 
                CGAL::canonicalize(tmp) 
                == 
                CGAL::canonicalize(d) );
}


template<class AT>
void test_univariate() {

    ::CGAL::set_pretty_mode(std::cout);

    typedef typename AT::Integer Integer;
   
    // testing univariate polynomials with integer coefficients
    typedef CGAL::Polynomial<Integer> int_Poly;
    typedef Integer NT;

    // special cases: 
    gcd_utcf_test(int_Poly(0),int_Poly(0),int_Poly(1));
    gcd_utcf_test(int_Poly(0),int_Poly(NT(4),NT(2)),int_Poly(NT(4),NT(2)));
    gcd_utcf_test(int_Poly(NT(4),NT(2)),int_Poly(0),int_Poly(NT(4),NT(2)));
  
    int_Poly d, f, g, h, result;
    int l;
    // lc is the first prime
    f = int_Poly(23,4,67109417);
    do{
        g = CGAL::CGALi::rand_Poly_int<Integer>(my_random.get_int(10,1000));
        result = CGAL::CGALi::gcd_utcf(f, g);
    }while(result.degree()!= 0);
    d = int_Poly(11, 2, 3, 1);
    gcd_utcf_test(f*d, g*d, d);

    // lc is the second prime
    f = int_Poly(23,4,67109431);
    do{
        g = CGAL::CGALi::rand_Poly_int<Integer>(my_random.get_int(10,1000));
        result = CGAL::CGALi::gcd_utcf(f, g);
    }while(result.degree()!= 0);
    d = int_Poly(11, 2, 3, 1);
    gcd_utcf_test(f*d, g*d, d);
  

    // random polynomials
    for(l=0;l<100;l++){
        do{
            f = CGAL::CGALi::rand_Poly_int<Integer>(my_random.get_int(10,1000));
            g = CGAL::CGALi::rand_Poly_int<Integer>(my_random.get_int(10,1000));
            result = CGAL::CGALi::gcd_utcf(f, g);     
        }
        while(result.degree()!= 0); 
   
        d = CGAL::CGALi::rand_Poly_int<Integer>(my_random.get_int(10,1000));
        int_Poly p1 = f*d;
        int_Poly p2 = g*d;
        gcd_utcf_test(p1, p2, d);
    }
  

  
// testing univariate polynomials with sqrt extension coefficients 
    typedef CGAL::Sqrt_extension<Integer  ,Integer> int_EXT_1;
    typedef CGAL::Sqrt_extension<int_EXT_1,Integer> int_EXT_2;  
    typedef CGAL::Polynomial<int_EXT_1> sqrt_Poly;

  
    sqrt_Poly e, m, n, result_sqrt;

    // special case, root is the first prime
    m = sqrt_Poly(int_EXT_1(Integer(232), Integer(2543), Integer(67109417)),
            int_EXT_1(Integer(34), Integer(25), Integer(67109417)), 
            int_EXT_1(Integer(4532124), Integer(2313), Integer(67109417)));
    n = sqrt_Poly(int_EXT_1(Integer(42), Integer(223), Integer(67109417)), 
            int_EXT_1(Integer(987144543), Integer(12125), Integer(67109417)));
  
    result_sqrt = CGAL::CGALi::gcd_utcf(m, n); 
    CGAL_postcondition_msg(result_sqrt.degree() == 0, 
            " wrong polynomials in sqrt special cases"); 
    e = sqrt_Poly(int_EXT_1(Integer(234), Integer(24341), Integer(67109417)), 
            int_EXT_1(Integer(42349), Integer(343251), Integer(67109417)));
  sqrt_gcd_utcf_test(m*e, n*e , e);


// special case, root is the second prime
    m = sqrt_Poly(int_EXT_1(Integer(2123323),Integer(23233),Integer(67109431)), 
            int_EXT_1(Integer(34), Integer(25), Integer(67109431)), 
            int_EXT_1(Integer(932345245), Integer(13323), Integer(67109431)));
    n = sqrt_Poly(int_EXT_1(Integer(434232), Integer(27823), Integer(67109431)),
            int_EXT_1(Integer(823178856), Integer(178825), Integer(67109431)));

    result_sqrt = CGAL::CGALi::gcd_utcf(m, n);  
    CGAL_postcondition_msg(result_sqrt.degree() == 0, 
            " wrong polynomials in sqrt special cases"); 
    e = sqrt_Poly(int_EXT_1(Integer(434214), Integer(21), Integer(67109431)), 
            int_EXT_1(Integer(9654349), Integer(22351), Integer(67109431)));
  sqrt_gcd_utcf_test(m*e, n*e , e);
 

    for(l=0;l<10;l++){
        do{
            m = CGAL::CGALi::rand_Poly_sqrt<int_EXT_1, Integer>
                (my_random.get_int(10,1000),Integer(789234));
            n = CGAL::CGALi::rand_Poly_sqrt<int_EXT_1, Integer>
                (my_random.get_int(10,1000),Integer(789234));
            result_sqrt = CGAL::CGALi::gcd_utcf(m, n);   
        }
        while(result_sqrt.degree()!= 0); 
   
        e = CGAL::CGALi::rand_Poly_sqrt<int_EXT_1, Integer>
            (my_random.get_int(10,1000),Integer(789234));
        sqrt_Poly q1 = m*e;
        sqrt_Poly q2 = n*e;
        sqrt_gcd_utcf_test(q1, q2, e);
    }	
}


int main(){
    
    // Set wrong rounding mode to test modular arithmetic 
    CGAL::Protect_FPU_rounding<true> pfr(CGAL_FE_UPWARD);

#ifdef CGAL_USE_LEDA
    test_univariate<CGAL::LEDA_arithmetic_kernel>();
#endif // CGAL_USE_LEDA    

#ifdef CGAL_USE_CORE  
    test_univariate<CGAL::CORE_arithmetic_kernel>();
#endif // Lis_HAVE_CORE
     
    return 0;
}


// EOF
