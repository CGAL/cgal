// Copyright (c) 2008 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
//
// Author(s)     : Michael Hemmer <hemmer@mpi-inf.mpg.de>


#ifndef CGAL_POLYNOMIAL_RESULTANT_H
#define CGAL_POLYNOMIAL_RESULTANT_H

// Modular arithmetic is slower, hence the default is 0
#ifndef CGAL_RESULTANT_USE_MODULAR_ARITHMETIC
#define CGAL_RESULTANT_USE_MODULAR_ARITHMETIC 0
#endif

#ifndef CGAL_RESULTANT_USE_DECOMPOSE
#define CGAL_RESULTANT_USE_DECOMPOSE 1
#endif 


#include <CGAL/basic.h>
#include <CGAL/Polynomial.h>

#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/Polynomial/Interpolator.h>
#include <CGAL/Polynomial/prs_resultant.h>
#include <CGAL/Polynomial/bezout_matrix.h>

#include <CGAL/Residue.h>
#include <CGAL/Modular_traits.h>
#include <CGAL/Chinese_remainder_traits.h>
#include <CGAL/primes.h>
#include <CGAL/Polynomial/Cached_extended_euclidean_algorithm.h>

namespace CGAL {


// The main function provided within this file is CGAL::internal::resultant(F,G),
// all other functions are used for dispatching. 
// The implementation uses interpolatation for multivariate polynomials
// Due to the recursive structuture of CGAL::Polynomial<Coeff> it is better 
// to write the function such that the inner most variabel is eliminated. 
// However,  CGAL::internal::resultant(F,G) eliminates the outer most variabel.
// This is due to backward compatibility issues with code base on EXACUS. 
// In turn CGAL::internal::resultant_(F,G) eliminates the innermost variable. 
 
// Dispatching
// CGAL::internal::resultant_decompose applies if Coeff is a Fraction
// CGAL::internal::resultant_modularize applies if Coeff is Modularizable
// CGAL::internal::resultant_interpolate applies for multivairate polynomials
// CGAL::internal::resultant_univariate selects the proper algorithm for IC 

// CGAL_RESULTANT_USE_DECOMPOSE ( default = 1 )
// CGAL_RESULTANT_USE_MODULAR_ARITHMETIC (default = 0 ) 

namespace internal{

template <class Coeff> 
inline Coeff resultant_interpolate( 
    const CGAL::Polynomial<Coeff>&, const CGAL::Polynomial<Coeff>& );
template <class Coeff> 
inline Coeff resultant_modularize( 
    const CGAL::Polynomial<Coeff>&, 
    const CGAL::Polynomial<Coeff>&, CGAL::Tag_true);
template <class Coeff> 
inline Coeff resultant_modularize( 
    const CGAL::Polynomial<Coeff>&, 
    const CGAL::Polynomial<Coeff>&, CGAL::Tag_false);
template <class Coeff> 
inline Coeff resultant_decompose( 
    const CGAL::Polynomial<Coeff>&, 
    const CGAL::Polynomial<Coeff>&, CGAL::Tag_true);
template <class Coeff> 
inline Coeff resultant_decompose( 
    const CGAL::Polynomial<Coeff>&, 
    const CGAL::Polynomial<Coeff>&, CGAL::Tag_false);
template <class Coeff> 
inline Coeff resultant_( 
    const CGAL::Polynomial<Coeff>&, const CGAL::Polynomial<Coeff>&);

template <class Coeff> 
inline Coeff resultant_univariate( 
    const CGAL::Polynomial<Coeff>& A, 
    const CGAL::Polynomial<Coeff>& B, 
    CGAL::Integral_domain_without_division_tag){ 
  return hybrid_bezout_subresultant(A,B,0);
}
template <class Coeff> 
inline Coeff resultant_univariate( 
    const CGAL::Polynomial<Coeff>& A, 
    const CGAL::Polynomial<Coeff>& B, CGAL::Integral_domain_tag){
  // this seems to help for for large polynomials 
  return prs_resultant_integral_domain(A,B);
}
template <class Coeff> 
inline Coeff resultant_univariate( 
    const CGAL::Polynomial<Coeff>& A, 
    const CGAL::Polynomial<Coeff>& B, CGAL::Unique_factorization_domain_tag){
  return prs_resultant_ufd(A,B);
}

template <class Coeff> 
inline Coeff resultant_univariate( 
    const CGAL::Polynomial<Coeff>& A, 
    const CGAL::Polynomial<Coeff>& B, CGAL::Field_tag){
  return prs_resultant_field(A,B);  
}

} // namespace internal

namespace internal{


template <class IC> 
inline IC 
resultant_interpolate( 
    const CGAL::Polynomial<IC>& F, 
    const CGAL::Polynomial<IC>& G){
  CGAL_precondition(CGAL::Polynomial_traits_d<CGAL::Polynomial<IC> >::d == 1);
    typedef CGAL::Algebraic_structure_traits<IC> AST_IC;
    typedef typename AST_IC::Algebraic_category Algebraic_category;
    return internal::resultant_univariate(F,G,Algebraic_category()); 
}

template <class Coeff_2> 
inline
CGAL::Polynomial<Coeff_2>  resultant_interpolate(
        const CGAL::Polynomial<CGAL::Polynomial<Coeff_2> >& F, 
        const CGAL::Polynomial<CGAL::Polynomial<Coeff_2> >& G){
    
    typedef CGAL::Polynomial<Coeff_2> Coeff_1;
    typedef CGAL::Polynomial<Coeff_1> POLY;
    typedef CGAL::Polynomial_traits_d<POLY> PT;
    typedef typename PT::Innermost_coefficient_type IC; 

    CGAL_precondition(PT::d >= 2);
    
    typename PT::Degree degree; 
    int maxdegree = degree(F,0)*degree(G,PT::d-1) + degree(F,PT::d-1)*degree(G,0); 

    typedef std::pair<IC,Coeff_2> Point; 
    std::vector<Point> points; // interpolation points  
    
   
    typename CGAL::Polynomial_traits_d<Coeff_1>::Degree  coeff_degree; 
    int i(-maxdegree/2);
    int deg_f(0);
    int deg_g(0);
    
   
    while((int) points.size() <= maxdegree + 1){
        i++;
        // timer1.start();
        Coeff_1 c_i(i);
        Coeff_1 Fat_i(typename PT::Evaluate()(F,c_i));
        Coeff_1 Gat_i(typename PT::Evaluate()(G,c_i));
        // timer1.stop();
        
        int deg_f_at_i = coeff_degree(Fat_i,0);
        int deg_g_at_i = coeff_degree(Gat_i,0);

        // std::cout << F << std::endl;
        // std::cout << Fat_i << std::endl;
        // std::cout << deg_f_at_i << " vs. " << deg_f << std::endl;
        if(deg_f_at_i >  deg_f ){
            points.clear();
            deg_f  = deg_f_at_i;
            CGAL_postcondition(points.size() == 0);
        } 

        if(deg_g_at_i >  deg_g){
            points.clear();
            deg_g  = deg_g_at_i;
            CGAL_postcondition(points.size() == 0);
        }
        
        if(deg_f_at_i ==  deg_f && deg_g_at_i ==  deg_g){
            // timer2.start();
            Coeff_2 res_at_i = resultant_interpolate(Fat_i, Gat_i);
            // timer2.stop();
            points.push_back(Point(IC(i),res_at_i));
            
            // std::cout << typename Polynomial_traits_d<Coeff_2>::Degree()(res_at_i) << std::endl ; 
        }      
    }
   
    // timer3.start();
    CGAL::internal::Interpolator<Coeff_1> interpolator(points.begin(),points.end());
    Coeff_1 result = interpolator.get_interpolant();
    // timer3.stop();

#ifndef CGAL_NDEBUG
    while((int) points.size() <= maxdegree + 3){
        i++;        

        Coeff_1 c_i(i);
        Coeff_1 Fat_i(typename PT::Evaluate()(F,c_i));
        Coeff_1 Gat_i(typename PT::Evaluate()(G,c_i));

        CGAL_assertion(coeff_degree(Fat_i,0) <= deg_f);
        CGAL_assertion(coeff_degree(Gat_i,0) <= deg_g);
        
        if(coeff_degree( Fat_i , 0) ==  deg_f && coeff_degree( Gat_i , 0 ) ==  deg_g){
            Coeff_2 res_at_i = resultant_interpolate(Fat_i, Gat_i);
            points.push_back(Point(IC(i), res_at_i));
        }
    }
    CGAL::internal::Interpolator<Coeff_1> 
      interpolator_(points.begin(),points.end());
    Coeff_1 result_= interpolator_.get_interpolant();
    
     // the interpolate polynomial has to be stable !
    CGAL_assertion(result_ == result); 
#endif 
    return result; 
}

template <class Coeff> 
inline
Coeff resultant_modularize( 
        const CGAL::Polynomial<Coeff>& F, 
        const CGAL::Polynomial<Coeff>& G, 
        CGAL::Tag_false){
    return resultant_interpolate(F,G);
}

template <class Coeff> 
inline
Coeff resultant_modularize( 
        const CGAL::Polynomial<Coeff>& F, 
        const CGAL::Polynomial<Coeff>& G, 
        CGAL::Tag_true){
    
  // Enforce IEEE double precision and to nearest before using modular arithmetic
  CGAL::Protect_FPU_rounding<true> pfr(CGAL_FE_TONEAREST);
  

    typedef Polynomial_traits_d<CGAL::Polynomial<Coeff> > PT;
    typedef typename PT::Polynomial_d Polynomial;
    
    typedef Chinese_remainder_traits<Coeff> CRT;
    typedef typename CRT::Scalar_type Scalar;


    typedef typename CGAL::Modular_traits<Polynomial>::Residue_type MPolynomial; 
    typedef typename CGAL::Modular_traits<Coeff>::Residue_type      MCoeff; 
        
    typename CRT::Chinese_remainder chinese_remainder; 
    typename CGAL::Modular_traits<Coeff>::Modular_image_representative inv_map;


    typename PT::Degree_vector                                     degree_vector; 
    typename CGAL::Polynomial_traits_d<MPolynomial>::Degree_vector mdegree_vector;

    bool solved = false; 
    int prime_index = 0; 
    int n = 0;
    Scalar p,q,pq,s,t; 
    Coeff R, R_old; 
    
    // CGAL::Timer timer_evaluate, timer_resultant, timer_cr; 
    
    do{
        MPolynomial mF, mG;
        MCoeff mR;
        //timer_evaluate.start();
        do{
            // select a prime number
            int current_prime = -1;
            prime_index++;
            if(prime_index >= 2000){
                std::cerr<<"primes in the array exhausted"<<std::endl;
                CGAL_assertion(false);
                current_prime = internal::get_next_lower_prime(current_prime);
            } else{
                current_prime = internal::primes[prime_index];
            }
            CGAL::Residue::set_current_prime(current_prime);
            
            mF = CGAL::modular_image(F);
            mG = CGAL::modular_image(G);
            
        }while( degree_vector(F) != mdegree_vector(mF) || 
                degree_vector(G) != mdegree_vector(mG));
        //timer_evaluate.stop();
        
        //timer_resultant.start();
        n++;
        mR = resultant_interpolate(mF,mG);
        //timer_resultant.stop();
        //timer_cr.start();
        if(n == 1){ 
            // init chinese remainder
            q =  CGAL::Residue::get_current_prime(); // implicit ! 
            R = inv_map(mR);
        }else{
            // continue chinese remainder
            p = CGAL::Residue::get_current_prime(); // implicit!  
            R_old  = R ;
//            chinese_remainder(q,Gs ,p,inv_map(mG_),pq,Gs);             
//            cached_extended_euclidean_algorithm(q,p,s,t);
            internal::Cached_extended_euclidean_algorithm
                <typename CRT::Scalar_type> ceea;
            ceea(q,p,s,t);
            pq =p*q;
            chinese_remainder(q,p,pq,s,t,R_old,inv_map(mR),R);
            q=pq;
        }
        solved = (R==R_old);
        //timer_cr.stop();       
    } while(!solved);
        
    //std::cout << "Time Evaluate   : " << timer_evaluate.time() << std::endl; 
    //std::cout << "Time Resultant  : " << timer_resultant.time() << std::endl; 
    //std::cout << "Time Chinese R  : " << timer_cr.time() << std::endl; 
    // CGAL_postcondition(R == resultant_interpolate(F,G));
    return R;
    // return resultant_interpolate(F,G);
}


template <class Coeff> 
inline
Coeff resultant_decompose( 
    const CGAL::Polynomial<Coeff>& F,
    const CGAL::Polynomial<Coeff>& G, 
    CGAL::Tag_false){
#if CGAL_RESULTANT_USE_MODULAR_ARITHMETIC
  typedef CGAL::Polynomial<Coeff> Polynomial; 
  typedef typename Modular_traits<Polynomial>::Is_modularizable Is_modularizable; 
  return resultant_modularize(F,G,Is_modularizable());
#else
  return  resultant_modularize(F,G,CGAL::Tag_false());
#endif
}

template <class Coeff> 
inline
Coeff resultant_decompose( 
        const CGAL::Polynomial<Coeff>& F, 
        const CGAL::Polynomial<Coeff>& G, 
        CGAL::Tag_true){  
    
    typedef Polynomial<Coeff> POLY;
    typedef typename Fraction_traits<POLY>::Numerator_type Numerator;
    typedef typename Fraction_traits<POLY>::Denominator_type Denominator;
    typename Fraction_traits<POLY>::Decompose decompose;
    typedef typename Numerator::NT RES;
    
    Denominator a, b;
    // F.simplify_coefficients(); not const 
    // G.simplify_coefficients(); not const 
    Numerator F0; decompose(F,F0,a);
    Numerator G0; decompose(G,G0,b);
    Denominator c = CGAL::ipower(a, G.degree()) * CGAL::ipower(b, F.degree());

    RES res0 =  CGAL::internal::resultant_(F0, G0);
    typename Fraction_traits<Coeff>::Compose comp_frac;
    Coeff res = comp_frac(res0, c);
    typename Algebraic_structure_traits<Coeff>::Simplify simplify;
    simplify(res);
    return res;
}


template <class Coeff> 
inline
Coeff resultant_( 
        const CGAL::Polynomial<Coeff>& F, 
        const CGAL::Polynomial<Coeff>& G){
#if CGAL_RESULTANT_USE_DECOMPOSE
    typedef CGAL::Fraction_traits<Polynomial<Coeff > > FT;
    typedef typename FT::Is_fraction Is_fraction; 
    return resultant_decompose(F,G,Is_fraction());
#else
    return resultant_decompose(F,G,CGAL::Tag_false());
#endif
}



template <class Coeff> 
inline
Coeff  resultant( 
        const CGAL::Polynomial<Coeff>& F_, 
        const CGAL::Polynomial<Coeff>& G_){
  // make the variable to be elimnated the innermost one.
    typedef CGAL::Polynomial_traits_d<CGAL::Polynomial<Coeff> > PT;
    CGAL::Polynomial<Coeff> F = typename PT::Move()(F_, PT::d-1, 0);
    CGAL::Polynomial<Coeff> G = typename PT::Move()(G_, PT::d-1, 0);
    return internal::resultant_(F,G);
}

} // namespace internal    
} //namespace CGAL



#endif // CGAL_POLYNOMIAL_RESULTANT_H
