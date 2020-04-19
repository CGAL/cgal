// Copyright (c) 2002-2008 Max-Planck-Institute Saarbruecken (Germany)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Dominik Huelse <dominik.huelse@gmx.de>
//                 Michael Hemmer <hemmer@mpi-inf.mpg.de>
//
// ============================================================================

// This is experimental code !

/*! \file CGAL/Polynomial/modular_gcd_utcf_with_wang.h
  provides gcd for Polynomials, based on Modular arithmetic
  with wang's algorithm as fallback.
*/


#ifndef CGAL_POLYNOMIAL_MODULAR_GCD_UTCF_WITH_WANG_H
#define CGAL_POLYNOMIAL_MODULAR_GCD_UTCF_WITH_WANG_H 1

#include <CGAL/basic.h>
#include <CGAL/Residue.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Scalar_factor_traits.h>
#include <CGAL/Chinese_remainder_traits.h>
#include <CGAL/Polynomial/modular_gcd_utils.h>
#include <CGAL/Polynomial/Cached_extended_euclidean_algorithm.h>
#include <CGAL/Timer.h>
#include <CGAL/Cache.h>
#include <CGAL/Polynomial/Wang_traits.h>

namespace CGAL {
namespace internal{

template <class NT> Polynomial<NT>
gcd_utcf_Integral_domain(Polynomial<NT>,Polynomial<NT>);


template <class NT>
Polynomial< Polynomial<NT> > modular_gcd_utcf_with_wang(
        const Polynomial< Polynomial<NT> >& FF1 ,
        const Polynomial< Polynomial<NT> >& FF2 ){
#ifndef NDEBUG
    std::cerr << "WARNING: still using internal::gcd_utcf_Integral_domain  "
              << std::endl;
    std::cerr << "TODO: fix mulitvariate modular_gcd_with_wang  "
              << std::endl;
#endif
    return gcd_utcf_Integral_domain(FF1, FF2);
}


// TODO: ALGORITHM M with Wang
template <class NT>
Polynomial<NT> modular_gcd_utcf_with_wang(
        const Polynomial<NT>& FF1_ ,
        const Polynomial<NT>& FF2_ ){

  // Enforce IEEE double precision and to nearest before using modular arithmetic
  CGAL::Protect_FPU_rounding<true> pfr(CGAL_FE_TONEAREST);

//    std::cout << "start modular_gcd_utcf_with_wang " << std::endl;

#ifdef CGAL_MODULAR_GCD_TIMER
    timer_init.start();
#endif


    typedef Polynomial<NT> Poly;
    typedef Polynomial_traits_d<Poly> PT;
    typedef CGAL::internal::Wang_traits<Poly> WT_POLY;


    typename WT_POLY::Wang wang;

    typedef typename PT::Innermost_coefficient_type IC;

    typename CGAL::Coercion_traits<Poly,IC>::Cast ictp;
    typename PT::Innermost_leading_coefficient ilcoeff;

    typedef Algebraic_extension_traits<IC> ANT;
    typename ANT::Denominator_for_algebraic_integers dfai;
    typename ANT::Normalization_factor nfac;

    // will play the role of content
    typedef typename CGAL::Scalar_factor_traits<Poly>::Scalar  Scalar;
    typedef typename CGAL::Modular_traits<Poly>::Residue_type   MPoly;
    typedef typename CGAL::Modular_traits<Scalar>::Residue_type MScalar;

    typedef Chinese_remainder_traits<Poly> CRT;
    typename CRT::Chinese_remainder chinese_remainder;

    Poly FF1 = FF1_;
    Poly FF2 = FF2_;

    if(FF1.is_zero()){
        if(FF2.is_zero()){
            //        std::cout<<"\nboth zero"<<std::endl;
            return Poly(1);
        }
        else{
            return CGAL::canonicalize(FF2);
        }
    }
    if(FF2.is_zero()){
        return CGAL::canonicalize(FF1);
    }
    if(FF1.degree() == 0 || FF2.degree() == 0){
        //      std::cout<<"\nconst polynomial"<<std::endl;
        return Poly(1);
    }

    // do we need this in case of wang?
    Poly F1 = CGAL::canonicalize(FF1);
    Poly F2 = CGAL::canonicalize(FF2);

    // in case IC is an algebraic extension it may happen, that
    // Fx=G*Hx is not possible if the coefficients are algebraic integers
    Poly tmp = F1+F2;
    IC denom = dfai(tmp.begin(),tmp.end());
    denom *= nfac(denom);

    Scalar denominator = scalar_factor(denom);

    // we use g_, since we hope that wang is not needed in this case
    Scalar f1 = scalar_factor(F1.lcoeff());  // ilcoeff(F1)
    Scalar f2 = scalar_factor(F2.lcoeff());  // ilcoeff(F2)
    Scalar g_ = scalar_factor(f1,f2);

    bool solved = false;
    int prime_index = -1;
    int n = 0; // number of lucky primes

    MScalar mg_;
    MPoly   mF1,mF2,mG_,mQ,mR;

    typename CRT::Scalar_type p,q,pq,s,t,m;
    Poly Gs,H1s,H2s, Gs_old;
    Poly Gsw(0), Gsw_old(0);
    bool wang_bool = false;

    typedef std::vector<int> Vector;
    Vector prs_degrees_old, prs_degrees_new;

    Poly result;
    int current_prime = -1;

    CGAL::Timer cr_timer, wang_timer;

#ifdef CGAL_MODULAR_GCD_TIMER
    timer_init.stop();
#endif


    while(!solved){
        do{
            //---------------------------------------
            //choose prime not deviding f1 or f2
            MScalar tmp1, tmp2;
            do{
                prime_index++;
                if(prime_index >= 2000){
                    std::cerr<<"primes exhausted"<<std::endl;
                    current_prime =
                        internal::get_next_lower_prime(current_prime);
                }
                else{
                    current_prime = internal::primes[prime_index];
                }
                CGAL_assertion(current_prime != -1);
                CGAL::Residue::set_current_prime(current_prime);
#ifdef CGAL_MODULAR_GCD_TIMER
                timer_image.start();
#endif
                tmp1 = modular_image(f1);
                tmp2 = modular_image(f2);
#ifdef CGAL_MODULAR_GCD_TIMER
                timer_image.stop();
#endif
            }
            // until prime does not devide both leading coeffs and denominator
            while(!(( tmp1 != 0 )
                            && ( tmp2 != 0 )
                            && (denominator%current_prime) != 0 ));

            // --------------------------------------
            // invoke gcd for current prime
#ifdef CGAL_MODULAR_GCD_TIMER
            timer_image.start();
#endif
            mg_ = CGAL::modular_image(g_);
            mF1 = CGAL::modular_image(FF1);
            mF2 = CGAL::modular_image(FF2);
#ifdef CGAL_MODULAR_GCD_TIMER
            timer_image.stop();
            // replace mG_ = CGAL::gcd (mF1,mF2)*MPoly(mg_); for multivariat
            timer_gcd.start();
#endif

            // compute gcd over Field[x]
            prs_degrees_new.clear();
            CGAL_precondition(mF1 != MPoly(0));
            CGAL_precondition(mF2 != MPoly(0));

            while (!mF2.is_zero()) {
                euclidean_division_obstinate(mF1, mF2, mQ, mR);
//                MPoly::euclidean_division(mF1, mF2, mQ, mR);
                mF1 = mF2; mF2 = mR;
                prs_degrees_new.push_back(mR.degree());
            }

            mF1 /= mF1.lcoeff();
            mG_ = mF1;

            if(n==0)
                prs_degrees_old = prs_degrees_new;

            mG_ = mG_ * MPoly(mg_);

#ifdef CGAL_MODULAR_GCD_TIMER
            timer_gcd.stop();
#endif

            //---------------------------------------
            // return if G is constant
            if (mG_ == MPoly(1)) return Poly(1);

            // use ordinary algorithm if prs sequence is too short
            // --------------------------------------
        }
        // repeat until mG_ degree is less equal the known bound
        // this is now tested by the prs degree sequence
        while(prs_degrees_old > prs_degrees_new);
        // check that everything went fine

        if( prs_degrees_old < prs_degrees_new ){
            if( n != 0 ) std::cout << "UNLUCKY PRIME !!"<< std::endl;

            // restart chinese remainder
            // ignore previous unlucky primes
            n=1;
            prs_degrees_old = prs_degrees_new;
        }else{
            CGAL_postcondition( prs_degrees_old == prs_degrees_new);
            n++; // increase number of lucky primes
        }

        typename CGAL::Modular_traits<Poly>::Modular_image_representative inv_map;


// ----------------------------- Chinese Remainder ---------------------
        cr_timer.start();
        if(n == 1){
            // init chinese remainder
            q =  CGAL::Residue::get_current_prime(); // implicit !
            Gs  = inv_map(mG_);
        }else{
            // continue chinese remainder
            p = CGAL::Residue::get_current_prime(); // implicit!
            Gs_old  = Gs ;
#ifdef CGAL_MODULAR_GCD_TIMER
            timer_CR.start();
#endif
            internal::Cached_extended_euclidean_algorithm < Scalar,3 > ceea;
            ceea(q,p,s,t);
            pq =p*q;
            chinese_remainder(q,p,pq,s,t,Gs,inv_map(mG_),Gs);
#ifdef CGAL_MODULAR_GCD_TIMER
            timer_CR.stop();
#endif
            q=pq;
        }

        try{
            // to catch error in case the extension is not correct yet.
            if( n != 1 && Gs == Gs_old ){
//                Poly r1,r2; NT dummy;
#ifdef CGAL_MODULAR_GCD_TIMER
                timer_division.start();
#endif
//                Gs = CGAL::canonicalize(Gs);
                typedef CGAL::Algebraic_structure_traits< Poly > ASTE_Poly;
                typename ASTE_Poly::Divides divides;


                FF1*=ictp(ilcoeff(Gs)*denom);
                FF2*=ictp(ilcoeff(Gs)*denom);

                bool div1=divides(Gs,FF1,H1s);
                bool div2=divides(Gs,FF2,H2s);

                if (div1 && div2){
                    solved = true;
                    result = Gs;
//                    std::cout << "number of primes used : "<< n << std::endl;
                }

                // this is the old code
//                Poly::euclidean_division(FF1,Gs,H1s,r1);  //TODO Divides
//                Poly::euclidean_division(FF2,Gs,H2s,r2);  //TODO Divides
//                if (r1.is_zero() && r2.is_zero()){
//                  solved = true;
//                  result = Gs;
//                  }



#ifdef CGAL_MODULAR_GCD_TIMER
                timer_division.stop();
#endif
            }
        }
        catch(...){}
        cr_timer.stop();


// ---------------------- wang ------------------------------
        // std::cout << cr_timer.time()<< " " << wang_timer.time()<< std::endl;
        if(!solved && cr_timer.time() > 50 * wang_timer.time() ){
            wang_timer.reset();
            cr_timer.reset();
            wang_timer.start();
            Gsw_old = Gsw;
            Poly   Gsw_;
            Scalar dummy;
#ifdef CGAL_MODULAR_GCD_TIMER
            timer_wang.start();
#endif
            wang_bool = wang(Gs,q,Gsw_,dummy);
#ifdef CGAL_MODULAR_GCD_TIMER
            timer_wang.stop();
#endif
            if(wang_bool){
                Gsw = Gsw_;
                try{ // if Gsw from wang is stable
                    if( Gsw == Gsw_old ) {
#ifdef CGAL_MODULAR_GCD_TIMER
                        timer_division.start();
#endif
//                        Poly r1,r2;
                        typedef CGAL::Algebraic_structure_traits
                            < Poly > ASTE_Poly;
                        typename ASTE_Poly::Divides divides;


                        FF1*=ictp(ilcoeff(Gsw)*denom);
                        FF2*=ictp(ilcoeff(Gsw)*denom);
                        bool div1=divides(Gsw,FF1,H1s);
                        bool div2=divides(Gsw,FF2,H2s);
                        // old code
//                        Poly::euclidean_division(FF1,Gsw,H1s,r1); //TODO Divides
//                        Poly::euclidean_division(FF2,Gsw,H2s,r2); //TODO Divides
#ifdef CGAL_MODULAR_GCD_TIMER
                        timer_division.stop();
#endif

                        if (div1 && div2){
                            solved = true;
                            result = Gsw;
                            //  std::cout << "number of primes used : "<< n << std::endl;
                        }

                        // old code
//                        if (r1.is_zero() && r2.is_zero()){
//                          solved = true;
//                          result = Gsw;
//                          }

                    }
                }catch(...){}


            }
            wang_timer.stop();
        }
    }
    return CGAL::canonicalize(result);

}
//#endif

}///namespace internal
}///namespace CGAL

#endif //#ifndef CGAL_POLYNOMIAL_MODULAR_GCD_UTCF_WITH_WANG_H

