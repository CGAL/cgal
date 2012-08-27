// Copyright (c) 2002-2008 Max-Planck-Institute Saarbruecken (Germany)
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
//
//
// Author(s)     :  Dominik Huelse <dominik.huelse@gmx.de>
//                  Michael Hemmer <mhemmer@uni-mainz.de>
//                 
// ============================================================================

/*! \file CGAL/Polynomial/modular_gcd_utcf_algorithm_M.h
  provides gcd for Polynomials, based on Modular arithmetic. 
*/


#ifndef CGAL_POLYNOMIAL_MODULAR_GCD_UTCF_ALGORITHM_M_H
#define CGAL_POLYNOMIAL_MODULAR_GCD_UTCF_ALGORITHM_M_H 1

#include <CGAL/basic.h>
#include <CGAL/Residue.h>
#include <CGAL/Polynomial/modular_gcd.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial/Cached_extended_euclidean_algorithm.h>
#include <CGAL/Scalar_factor_traits.h>
#include <CGAL/Chinese_remainder_traits.h>
#include <CGAL/Cache.h>
#include <CGAL/Real_timer.h>

// algorithm M for integer polynomials, without denominator bound

namespace CGAL {

namespace internal{
template <class NT> Polynomial<NT> gcd_utcf_UFD(Polynomial<NT>,Polynomial<NT>);


template <class NT> 
Polynomial< Polynomial<NT> > modular_gcd_utcf_algorithm_M(
        const Polynomial< Polynomial<NT> >& FF1 ,
        const Polynomial< Polynomial<NT> >& FF2 ){
    return gcd_utcf_UFD(FF1, FF2);
}

template <class NT> 
Polynomial<NT> modular_gcd_utcf_algorithm_M(
        const Polynomial<NT>& FF1 ,
        const Polynomial<NT>& FF2 ){

  // Enforce IEEE double precision and to nearest before using modular arithmetic
  CGAL::Protect_FPU_rounding<true> pfr(CGAL_FE_TONEAREST);

//    std::cout << "start modular_gcd_utcf_algorithm_M " << std::endl;
#ifdef CGAL_MODULAR_GCD_TIMER
    timer_init.start();
#endif
    typedef Polynomial<NT> Poly;

    // will paly the role of content
    typedef typename CGAL::Scalar_factor_traits<Poly>::Scalar  Scalar;

    typedef typename CGAL::Modular_traits<Poly>::Residue_type   MPoly;
    typedef typename CGAL::Modular_traits<Scalar>::Residue_type MScalar;
    
    typedef Chinese_remainder_traits<Poly> CRT;
    typename CRT::Chinese_remainder chinese_remainder; 
    
    CGAL::Real_timer timer; 
    
 
    if(FF1.is_zero()){
        if(FF2.is_zero()){
            return Poly(1);// TODO: return 0 for CGAL 
        }
        else{
            //      std::cout<<"\nFF1 is zero"<<std::endl;

            return CGAL::canonicalize(FF2);
        }
    }
    if(FF2.is_zero()){
        return CGAL::canonicalize(FF1);
    }
    if(FF1.degree() == 0 || FF2.degree() == 0){
        Poly result;
        result = Poly(CGAL::gcd(FF1.content(),FF2.content()));
        return CGAL::canonicalize(result);
    }
        
    Poly F1 = CGAL::canonicalize(FF1);
    Poly F2 = CGAL::canonicalize(FF2);
       
    Scalar f1 = scalar_factor(F1.lcoeff());  // ilcoeff(F1) 
    Scalar f2 = scalar_factor(F2.lcoeff());  // ilcoeff(F2) 
    Scalar g_ = scalar_factor(f1,f2);
    
    //std::cout <<" g_   : "<< g_ << std::endl;
    
    bool solved = false;
    int prime_index = -1;
    int n = 0; // number of lucky primes 
    int degree_F1 = F1.degree();
    int degree_F2 = F2.degree();
    int degree_e = (std::min)(degree_F1,degree_F2);

    MScalar mg_;
    MPoly   mF1,mF2,mG_;

    typename CRT::Scalar_type p,q,pq,s,t;
    Poly Gs,H1s,H2s, Gs_old; // s =^ star 
#ifdef CGAL_MODULAR_GCD_TIMER
    timer_init.stop();
#endif

    while(!solved){
        do{
            //---------------------------------------
            //choose prime not deviding f1 or f2
            MScalar tmp1, tmp2;
            do{
                int current_prime = -1;
                prime_index++;
		if(prime_index >= 2000){
                    std::cerr<<"primes in the array exhausted"<<std::endl;
                    current_prime = internal::get_next_lower_prime(current_prime);
                }
                else{
                    current_prime = internal::primes[prime_index];
                }
                CGAL::Residue::set_current_prime(current_prime);
#ifdef CGAL_MODULAR_GCD_TIMER
                timer_image.start();
#endif
                tmp1 = CGAL::modular_image(f1);
                tmp2 = CGAL::modular_image(f2);
#ifdef CGAL_MODULAR_GCD_TIMER
                timer_image.stop();
#endif
            }
            while(!(( tmp1 != 0 ) && ( tmp2 != 0 )));
            // --------------------------------------
            // invoke gcd for current prime
#ifdef CGAL_MODULAR_GCD_TIMER
            timer_image.start();
#endif
            mg_ = CGAL::modular_image(g_);
            mF1 = CGAL::modular_image(F1);
            mF2 = CGAL::modular_image(F2);
#ifdef CGAL_MODULAR_GCD_TIMER
            timer_image.stop();
            // replace mG_ = CGAL::gcd (mF1,mF2)*MPoly(mg_); for multivariat
            timer_gcd.start();
#endif
            mG_ = CGAL::gcd(mF1,mF2)*MPoly(mg_);
#ifdef CGAL_MODULAR_GCD_TIMER
            timer_gcd.stop();
#endif
            
            //mH1 = CGAL::integral_div(mF1,mG_);
            //mH2 = CGAL::integral_div(mF2,mG_);
            //---------------------------------------
            // return if G is constant 
            if (mG_ == MPoly(1)) return Poly(1);
            // --------------------------------------
        }// repeat until mG_ degree is less equal the known bound
         // check prime 
        while( mG_.degree() > degree_e);
       
        if( mG_.degree() < degree_e ){
            if( n != 0 ) std::cout << "UNLUCKY PRIME !!"<< std::endl; 

            // restart chinese remainder 
            // ignore previous unlucky primes
            n=1; 
            degree_e= mG_.degree();
        }else{
            CGAL_postcondition( mG_.degree() == degree_e);
            n++; // increase number of lucky primes
        }
     
        // --------------------------------------
        // try chinese remainder
        
//        std::cout <<" chinese remainder round :" << n << std::endl; 
        typename CGAL::Modular_traits<Poly>::Modular_image_representative inv_map;
        if(n == 1){ 
            // init chinese remainder
            q =  CGAL::Residue::get_current_prime(); // implicit ! 
            Gs_old  = Gs  = inv_map(mG_);
            
            //H1s_old = H1s = inv_map(mH1);
            //H2s_old = H2s = inv_map(mH2);
        }else{
            // continue chinese remainder
            
            p = CGAL::Residue::get_current_prime(); // implicit! 
             
            Gs_old  = Gs ;
            //H1s_old = H1s ;
            //H2s_old = H2s ;
#ifdef CGAL_MODULAR_GCD_TIMER
            timer_CR.start();
#endif
//            chinese_remainder(q,Gs ,p,inv_map(mG_),pq,Gs);             
//            cached_extended_euclidean_algorithm(q,p,s,t);
            internal::Cached_extended_euclidean_algorithm
              <typename CRT::Scalar_type, 1> ceea;
            ceea(q,p,s,t);
            pq =p*q;
            chinese_remainder(q,p,pq,s,t,Gs,inv_map(mG_),Gs);
#ifdef CGAL_MODULAR_GCD_TIMER
            timer_CR.stop();
#endif
            q=pq;
        }

        try{
            if( n != 1 && Gs == Gs_old ){
                Poly r1,r2; 
#ifdef CGAL_MODULAR_GCD_TIMER
                timer_division.start();
#endif 

                typedef CGAL::Algebraic_structure_traits< Poly > ASTE_Poly;
                typename ASTE_Poly::Divides divides;

                bool div1=divides(Gs,g_*F1,H1s);
                bool div2=divides(Gs,g_*F2,H2s);
                if (div1 && div2){
                    solved = true; 
                }
                // this is the old code
//                 NT dummy; 
//                Poly::euclidean_division(g_*F1,Gs,H1s,r1);
//                Poly::euclidean_division(g_*F2,Gs,H2s,r2);
//                if (r1.is_zero() && r2.is_zero())
//                      solved = true;       
      
#ifdef CGAL_MODULAR_GCD_TIMER
                timer_division.stop();
#endif
//                std::cout << "number of primes used : "<< n << std::endl;
            } // end while
           
        }catch(...){}
 
    }
    
    
    //TODO CGAL: change this to multivariat content
//    Scalar scalar_content_f1 = scalar_factor(FF1);
//    Scalar scalar_content_f2 = scalar_factor(FF2);
//    Scalar scalar_content_gcd = CGAL::gcd(scalar_content_f1,scalar_content_f2);
//    Poly result = CGAL::canonicalize(Gs)*Poly(scalar_content_gcd);
//    return result; 

    return CGAL::canonicalize(Gs);
    
}

}  // namespace internal

}  // namespace CGAL

#endif // CGAL_POLYNOMIAL_MODULAR_GCD_UTCF_ALGORITHM_M_H
