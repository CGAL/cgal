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

/*! \file CGAL/Polynomial/modular_gcd_utcf_dfai.h
  provides gcd for Polynomials, based on Modular arithmetic.
*/
#ifndef CGAL_POLYNOMIAL_MODULAR_GCD_UTCF_DFAI_H
#define CGAL_POLYNOMIAL_MODULAR_GCD_UTCF_DFAI_H 1

#include <CGAL/basic.h>
#include <CGAL/Scalar_factor_traits.h>
#include <CGAL/Chinese_remainder_traits.h>
#include <CGAL/Algebraic_extension_traits.h>
#include <CGAL/primes.h>
#include <CGAL/Residue.h>
#include <CGAL/Polynomial/modular_gcd.h>
#include <CGAL/Polynomial/Cached_extended_euclidean_algorithm.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial/modular_gcd_utils.h>
#include <CGAL/Cache.h>
#include <CGAL/Real_timer.h>



namespace CGAL {

namespace internal{

template <class NT> Polynomial<NT>
gcd_utcf_Integral_domain(Polynomial<NT>,Polynomial<NT>);


template <class NT>
Polynomial< Polynomial<NT> > modular_gcd_utcf_dfai(
        const Polynomial< Polynomial<NT> >& FF1 ,
        const Polynomial< Polynomial<NT> >& FF2 ){
    return internal::gcd_utcf_Integral_domain(FF1, FF2);
}

// algorithm just computes Gs using the denominator for
// algebraic integers  and checks it using pseudo division.
template <class NT>
Polynomial<NT> modular_gcd_utcf_dfai(
        const Polynomial<NT>& FF1_ ,
        const Polynomial<NT>& FF2_ ){

  // Enforce IEEE double precision and to nearest before using modular arithmetic
  CGAL::Protect_FPU_rounding<true> pfr(CGAL_FE_TONEAREST);

    typedef Polynomial<NT> Poly;
    Poly FF1 = FF1_;
    Poly FF2 = FF2_;

    typedef Polynomial_traits_d<Poly> PT;
    typedef typename PT::Innermost_coefficient_type IC;

    typename Coercion_traits<Poly,IC>::Cast ictp;
    typename PT::Construct_innermost_coefficient_const_iterator_range range;
    typename PT::Innermost_leading_coefficient ilcoeff;

    typedef  Algebraic_extension_traits<IC> ANT;
    typename ANT::Denominator_for_algebraic_integers dfai;
    typename ANT::Normalization_factor nfac;

    typedef typename CGAL::Scalar_factor_traits<Poly>  SFT;
    typedef typename SFT::Scalar Scalar;

    typedef typename CGAL::Modular_traits<Poly>::Residue_type   MPoly;
    typedef typename CGAL::Modular_traits<Scalar>::Residue_type MScalar;

    typedef Chinese_remainder_traits<Poly> CRT;
    typename CRT::Chinese_remainder chinese_remainder;

    CGAL::Real_timer timer;

    if(FF1.is_zero()){
        if(FF2.is_zero()){
            return Poly(1);
        }
        else{
            return CGAL::canonicalize(FF2);
        }
    }
    if(FF2.is_zero()){;
        return CGAL::canonicalize(FF1);
    }

    if(FF1.degree() == 0 || FF2.degree() == 0){
        return Poly(1);
    }

#ifdef CGAL_MODULAR_GCD_TIMER
    timer_init.start();
#endif

    Poly F1 = CGAL::canonicalize(FF1);
    Poly F2 = CGAL::canonicalize(FF2);


    // This is the most important part of the (dfai) algorithm, it computes the
    // multiplictive denominator bound according to the algebraic extension
    // This is needed to guarantee a termination of the Chenese Remainder,
    // i.e. ensures that Gs  can be expressed in terms of algebraic integers
    // dfai = denominator for algebraic integers

    IC denom;
    {
        Poly tmp = F1+F2;
        denom = dfai(range(tmp).first, range(tmp).second);
    }
    denom *= nfac(denom);

    Scalar denominator = scalar_factor(denom);
    //std::cout <<" F1*denom*nafc: " << F1 <<std::endl;
    //std::cout <<" F2*denom*nfac: " << F2 <<std::endl;

    Scalar f1 = scalar_factor(F1.lcoeff());  // ilcoeff(F1)
    Scalar f2 = scalar_factor(F2.lcoeff());  // ilcoeff(F2)
    Scalar g_ = scalar_factor(f1,f2) * denominator;

    bool solved = false;
    int prime_index = -1;

    int n = 0; // number of lucky primes

    MScalar mg_;
    MPoly   mF1,mF2,mG_,mQ,mR;

    typename CRT::Scalar_type p,q,pq,s,t;
    Poly Gs,H1s,H2s,Gs_old; // s =^ star

#ifdef CGAL_MODULAR_GCD_TIMER
    timer_init.stop();
#endif

    typedef std::vector<int> Vector;
    Vector prs_degrees_old, prs_degrees_new;
    int current_prime = -1;


    while(!solved){
        do{
            //---------------------------------------
            //choose prime not deviding f1 or f2
            MScalar tmp1, tmp2;
            do{
                prime_index++;
                if(prime_index >= 2000){
                    std::cerr<<"primes exhausted"<<std::endl;
//                    return CGAL::internal::gcd_utcf_Integral_domain(FF1,FF2);
                    current_prime = internal::get_next_lower_prime(current_prime);
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
#endif

#ifdef CGAL_MODULAR_GCD_TIMER
            timer_gcd.start();
#endif

            // compute gcd over Field[x]

            prs_degrees_new.clear();
            CGAL_precondition(mF1 != MPoly(0));
            CGAL_precondition(mF2 != MPoly(0));
            while (!mF2.is_zero()) {
//                MPoly::euclidean_division(mF1, mF2, mQ, mR);
                euclidean_division_obstinate(mF1, mF2, mQ, mR);
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

            // return if G is constant
            if (mG_ == MPoly(1)) return Poly(1);
            // use ordinary algorithm if prs sequence is too short
            //  if (prs_degrees_new.size() <= 2)
            //    return internal::gcd_utcf_Integral_domain(F1, F2);
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

        // --------------------------------------
        // try chinese remainder
        typename CGAL::Modular_traits<Poly>::Modular_image_representative inv_map;
        if(n == 1){
            // init chinese remainder
            q =  CGAL::Residue::get_current_prime(); // implicit !
            Gs_old  = Gs  = inv_map(mG_);
        }else{
            // continue chinese remainder
            p = CGAL::Residue::get_current_prime(); // implicit!
            Gs_old  = Gs ;
#ifdef CGAL_MODULAR_GCD_TIMER
            timer_CR.start();
#endif
            internal::Cached_extended_euclidean_algorithm <Scalar, 2> ceea;
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
                Poly r1,r2; NT dummy;
#ifdef CGAL_MODULAR_GCD_TIMER
                timer_division.start();
#endif
                typedef CGAL::Algebraic_structure_traits< Poly >
                    ASTE_Poly;
                typename ASTE_Poly::Divides divides;


                FF1*=ictp(ilcoeff(Gs)*denom);
                FF2*=ictp(ilcoeff(Gs)*denom);
                bool div1=divides(Gs,FF1,H1s);
                bool div2=divides(Gs,FF2,H2s);


                //This is the old code:
//                Poly::pseudo_division(F1,Gs,H1s,r1,dummy);
//                Poly::pseudo_division(F2,Gs,H2s,r2,dummy);
                if (div1 && div2){
                    solved = true;
//                    std::cout << "number of primes used : "<< n << std::endl;
                }

//                if (r1.is_zero() && r2.is_zero()){
//                    solved = true;
//                    std::cout << "number of primes used : "<< n << std::endl;
//                }

#ifdef CGAL_MODULAR_GCD_TIMER
                timer_division.stop();
#endif
            }
        }
        catch(...){}
    }

    return CGAL::canonicalize(Gs);

}

} // namespace internal
} // namespace CGAL

#endif // CGAL_POLYNOMIAL_MODULAR_GCD_UTCF_DFAI_H
