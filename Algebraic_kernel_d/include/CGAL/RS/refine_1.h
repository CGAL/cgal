// Copyright (c) 2006-2008 Inria Lorraine (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
// Author: Luis Peñaranda <luis.penaranda@loria.fr>

#ifndef CGAL_RS_REFINE_1_H
#define CGAL_RS_REFINE_1_H

#include <gmp.h>
#include <mpfr.h>
#include <mpfi.h>
#include <CGAL/RS/polynomial_1.h>
#include <CGAL/RS/algebraic_1.h>
#include <CGAL/RS/polynomial_1_utils.h>
#include <CGAL/RS/dyadic.h>
#include <CGAL/RS/sign_1.h>
#include <CGAL/assertions.h>

namespace CGAL{

class Algebraic_1;

// bisects a's isolation interval in point p, returns a
// Comparison_result, which is the comparison between a and p
// precondition: p belongs to a's isolation interval
template <class _Gcd_policy>
Comparison_result bisect_at(const Algebraic_1 &a,mpfr_srcptr p){
        typedef _Gcd_policy     Gcd;
        Sign sl,sp;
        int round;
        sp=RSSign::signat(sfpart_1<Gcd>()(a.pol()),p);
        if(sp==ZERO){
                // set a=[p,p]
                mpfr_set_prec(&(a.mpfi()->left),mpfr_get_prec(p));
                round=mpfr_set(&(a.mpfi()->left),p,GMP_RNDN);
                CGAL_assertion(!round);
                mpfr_set_prec(&(a.mpfi()->right),mpfr_get_prec(p));
                round=mpfr_set(&(a.mpfi()->right),p,GMP_RNDN);
                CGAL_assertion(!round);
                a.set_lefteval(ZERO);
                return EQUAL;
        }
        sl=a.lefteval();
        if(sp==sl){
                // set a=[p,a_right]
                mpfr_set_prec(&(a.mpfi()->left),mpfr_get_prec(p));
                round=mpfr_set(&(a.mpfi()->left),p,GMP_RNDN);
                CGAL_assertion(!round);
                return LARGER;
        }else{
                // set a=[a_left,p]
                mpfr_set_prec(&(a.mpfi()->right),mpfr_get_prec(p));
                round=mpfr_set(&(a.mpfi()->right),p,GMP_RNDN);
                CGAL_assertion(!round);
                return SMALLER;
        }
}

// refines b until having either only one point in common with a or
// all points inside a; the result will be zero when all points of
// b result to be in a, negative when a<b and positive when a>b
// precondition: a and b overlap
template <class _Gcd_policy>
int bisect_at_endpoints(const Algebraic_1 &a,Algebraic_1 &b){
        typedef _Gcd_policy     Gcd;
        CGAL_precondition(a.overlaps(b));
        if(mpfr_cmp(&(a.mpfi()->left),&(b.mpfi()->left))>0){
                Comparison_result refinement=
                        bisect_at<Gcd>(b,&(a.mpfi()->left));
                if(refinement==EQUAL)
                        return 0;
                if(refinement==SMALLER)
                        return 1;
        }
        if(mpfr_cmp(&(a.mpfi()->right),&(b.mpfi()->right))<0 &&
                        bisect_at<Gcd>(b,&(a.mpfi()->right))==LARGER)
                return -1;
        return 0;
}

// bisect n times an interval; this function returns the number of
// refinements made
template <class _Gcd_policy>
int bisect_n(const Algebraic_1 &a,unsigned long n=1){
        typedef _Gcd_policy     Gcd;
        Sign sl,sc;
        mp_prec_t pl,pc;
        mpfr_t center;
        unsigned long i;
        int round;

        sl=a.lefteval();
        if(sl==ZERO)
                return 0;
        pl=mpfr_get_prec(&(a.mpfi()->left));
        pc=mpfr_get_prec(&(a.mpfi()->right));
        pc=(pl>pc?pl:pc)+(mp_prec_t)n;
        mpfr_init2(center,pc);
        round=mpfr_prec_round((mpfr_ptr)&(a.mpfi()->left),pc,GMP_RNDN);
        CGAL_assertion(!round);
        round=mpfr_prec_round((mpfr_ptr)&(a.mpfi()->right),pc,GMP_RNDN);
        CGAL_assertion(!round);
        for(i=0;i<n;++i){
                round=mpfr_add(
                        center,
                        &(a.mpfi()->left),
                        &(a.mpfi()->right),
                        GMP_RNDN);
                CGAL_assertion(!round);
                round=mpfr_div_2ui(center,center,1,GMP_RNDN);
                CGAL_assertion(!round);
                sc=RSSign::signat(sfpart_1<Gcd>()(a.pol()),center);
                if(sc==ZERO){   // we have a root
                        round=mpfr_set((mpfr_ptr)&(a.mpfi()->left),
                                       center,
                                       GMP_RNDN);
                        CGAL_assertion(!round);
                        mpfr_swap((mpfr_ptr)&(a.mpfi()->right),center);
                        a.set_lefteval(ZERO);
                        break;
                }
                if(sc==sl)
                        mpfr_swap((mpfr_ptr)&(a.mpfi()->left),center);
                else
                        mpfr_swap((mpfr_ptr)&(a.mpfi()->right),center);
        }
        mpfr_clear(center);
        return i;
}

// The same bisection function, but with dyadic numbers. The difference
// is that dyadics will use the exact amount of bits needed, without
// allocating a big amount of memory. Note that the dyadic numbers are
// implemented as mpfrs, this implies that no conversion is needed.
template <class _Gcd_policy>
int bisect_n_dyadic(const Algebraic_1 &a,unsigned long n=1){

        typedef _Gcd_policy     Gcd;
        Sign sl,sc;
        CGALRS_dyadic_t center;
        unsigned  long i;
        sl=a.lefteval();
        if(sl==ZERO)
                return 0;
        CGALRS_dyadic_init(center);
        for(i=0;i<n;++i){
                CGALRS_dyadic_midpoint(center,a.left(),a.right());
                sc=RSSign::signat(sfpart_1<Gcd>()(a.pol()),center);
                if(sc==ZERO){   // we have a root
                        CGALRS_dyadic_set((CGALRS_dyadic_ptr)a.left(),center);
                        CGALRS_dyadic_swap((CGALRS_dyadic_ptr)a.right(),center);
                        a.set_lefteval(ZERO);
                        break;
                }
                if(sc==sl)
                        CGALRS_dyadic_swap((CGALRS_dyadic_ptr)a.left(),center);
                else
                        CGALRS_dyadic_swap((CGALRS_dyadic_ptr)a.right(),center);
        }
        CGALRS_dyadic_clear(center);
        return i;
}

// refine an interval, by bisecting it, until having a size smaller than 2^(-s)
template <class _Gcd_policy>
int bisect(const Algebraic_1 &a,mp_exp_t s){
        typedef _Gcd_policy     Gcd;
        long ed;
        mpfr_t d;       // interval size
        mpfr_init(d);
        mpfr_sub(d,a.right(),a.left(),GMP_RNDU);
        mpfr_log2(d,d,GMP_RNDU);
        ed=1+mpfr_get_si(d,GMP_RNDU);
        // we found ed such that the interval size is between 2^(ed-1) and 2^ed
        mpfr_clear(d);
        // following tests, dyadic bisection is a bit faster than mpfr one
        return bisect_n_dyadic<Gcd>(a,s-ed);
}

// the four following functions are used to apply quadratic interval
// refinement (Abbott, 2006)

// calculate the value of kappa
//--------------------------------------------------
// unsigned long calc_kappa(const RS_polynomial_1 &p,
//                 CGALRS_dyadic_ptr f_lo,CGALRS_dyadic_ptr f_hi,
//                 unsigned n){
//         unsigned long ret;
//         CGALRS_dyadic_t temp;
//         mpfr_t numerator,denominator;
//         mpfr_inits2(MPFR_PREC_MIN,numerator,denominator,NULL);
//         CGALRS_dyadic_init(temp);
//
//         CGALRS_dyadic_sub(temp,f_lo,f_hi);
//         CGALRS_dyadic_get_exactfr(denominator,temp);
//         CGALRS_dyadic_mul_2exp(temp,f_lo,n);
//         CGALRS_dyadic_get_exactfr(numerator,temp);
//
//         mpfr_div(numerator,numerator,denominator,GMP_RNDN);
//         ret=1+mpfr_get_ui(numerator,GMP_RNDN);
//
//         CGALRS_dyadic_clear(temp);
//         mpfr_clears(numerator,denominator,NULL);
//         return ret;
// }
//--------------------------------------------------
unsigned long calc_kappa(const RS_polynomial_1 &p,
                CGALRS_dyadic_ptr f_lo,CGALRS_dyadic_ptr f_hi,
                unsigned n){
        unsigned long ret;
        //--------------------------------------------------
        // CGALRS_dyadic_t numerator,denominator;
        // CGALRS_dyadic_init(numerator);
        // CGALRS_dyadic_init(denominator);
        //--------------------------------------------------

        CGALRS_dyadic_sub(f_hi,f_lo,f_hi);
        CGALRS_dyadic_mul_2exp(f_lo,f_lo,n);

        mpfr_div(f_lo,f_lo,f_hi,GMP_RNDN);
        ret=mpfr_get_ui(f_lo,GMP_RNDN);

        //--------------------------------------------------
        // CGALRS_dyadic_clear(numerator);
        // CGALRS_dyadic_clear(denominator);
        //--------------------------------------------------
        return ret;
}

// calculate kappa using doubles; faster but less accurate
//--------------------------------------------------
// inline unsigned calc_kappa(const RS_polynomial_1 &p,
//                 double f_lo,double f_hi,
//                 unsigned n){
//         return(1+(int)nearbyint(f_lo*pow(2,n)/(f_lo-f_hi)));
// }
//--------------------------------------------------

// returns 0 for success, -1 for failure and 1 if it finds the root;
// the "N" of the paper is represented here as 2^n
int refine_interval_by_factor(CGALRS_dyadic_ptr x_lo,CGALRS_dyadic_ptr x_hi,
                CGALRS_dyadic_ptr f_lo,CGALRS_dyadic_ptr f_hi,
                const RS_polynomial_1 &p,unsigned n){
        Sign sl=RSSign::signat(p,x_lo);
        unsigned long kappa;
        CGALRS_dyadic_t xkappa,xside,w;
        Sign s1,s2;
        kappa=calc_kappa(p,f_lo,f_hi,n);
        //kappa=calc_kappa(p,
        //                 CGALRS_dyadic_get_d(f_lo),
        //                 CGALRS_dyadic_get_d(f_hi),
        //                 n);
        if(n==2){       // N==4: bisect twice
                unsigned b[2],inter;
                // we will use first xkappa as bisection point
                CGALRS_dyadic_init(xkappa);
                for(int i=0;i<2;++i){
                        CGALRS_dyadic_midpoint(xkappa,x_lo,x_hi);
                        if((s1=RSSign::signat(p,xkappa))==ZERO){        //exact
                                CGALRS_dyadic_set(x_lo,xkappa);
                                CGALRS_dyadic_swap(x_hi,xkappa);
                                CGALRS_dyadic_clear(xkappa);
                                return 1;
                        }
                        if(sl==s1){     // we take the right half
                                b[i]=2-i;
                                CGALRS_dyadic_swap(x_lo,xkappa);
                        }else{  // we take the left half
                                b[i]=0;
                                CGALRS_dyadic_swap(x_hi,xkappa);
                        }
                }
                CGALRS_dyadic_clear(xkappa);
                inter=b[0]+b[1];
                switch(inter){
                        case 0:p.inexact_eval_mpfr(f_hi,x_hi);break;
                        case 3:p.inexact_eval_mpfr(f_lo,x_lo);break;
                        default:p.inexact_eval_mpfr(f_hi,x_hi);
                                p.inexact_eval_mpfr(f_lo,x_lo);
                                break;
                }
                if(inter+1==kappa)
                        return 0;       // success
                else
                        return -2;      // failure (but we refined anyway)
        }else{  // N!=4
                // calculate w, the width of every interval
                CGALRS_dyadic_init(w);
                CGALRS_dyadic_sub(w,x_hi,x_lo);
                CGALRS_dyadic_div_2exp(w,w,n);
                // calculate x[kappa]
                CGALRS_dyadic_init_set(xkappa,x_lo);
                CGALRS_dyadic_addmul_ui(xkappa,w,kappa);
                if((s1=RSSign::signat(p,xkappa))==ZERO){
                        CGALRS_dyadic_set(x_lo,xkappa);
                        CGALRS_dyadic_swap(x_hi,xkappa);
                        CGALRS_dyadic_clear(xkappa);
                        CGALRS_dyadic_clear(w);
                        return 1;       // exact root
                }
                if(s1==sl){
                        // calculate x[kappa+1]
                        CGALRS_dyadic_init(xside);
                        CGALRS_dyadic_add(xside,xkappa,w);
                        if((s2=RSSign::signat(p,xside))==ZERO){
                                CGALRS_dyadic_set(x_lo,xside);
                                CGALRS_dyadic_swap(x_hi,xside);
                                CGALRS_dyadic_clear(xkappa);
                                CGALRS_dyadic_clear(w);
                                CGALRS_dyadic_clear(xside);
                                return 1;       // exact root
                        }
                        if(s2==sl){
                                CGALRS_dyadic_clear(xkappa);
                                CGALRS_dyadic_clear(w);
                                CGALRS_dyadic_clear(xside);
                                return -1;      // failure
                        }else{  // s2!=sl && s2!=0
                                CGALRS_dyadic_set(x_lo,xkappa);
                                CGALRS_dyadic_swap(x_hi,xside);
                                p.inexact_eval_mpfr(f_lo,x_lo);
                                p.inexact_eval_mpfr(f_hi,x_hi);
                                CGALRS_dyadic_clear(xkappa);
                                CGALRS_dyadic_clear(w);
                                CGALRS_dyadic_clear(xside);
                                return 0;       // success
                        }
                }else{  // s1!=sl && s1!=ZERO
                        // we calculate x[kappa-1];
                        CGALRS_dyadic_init(xside);
                        CGALRS_dyadic_sub(xside,xkappa,w);
                        if((s2=RSSign::signat(p,xside))==ZERO){
                                CGALRS_dyadic_set(x_lo,xside);
                                CGALRS_dyadic_swap(x_hi,xside);
                                CGALRS_dyadic_clear(xkappa);
                                CGALRS_dyadic_clear(w);
                                CGALRS_dyadic_clear(xside);
                                return 1;       // exact root
                        }
                        if(s2==sl){
                                CGALRS_dyadic_swap(x_lo,xside);
                                CGALRS_dyadic_swap(x_hi,xkappa);
                                p.inexact_eval_mpfr(f_lo,x_lo);
                                p.inexact_eval_mpfr(f_hi,x_hi);
                                CGALRS_dyadic_clear(xkappa);
                                CGALRS_dyadic_clear(w);
                                CGALRS_dyadic_clear(xside);
                                return 0;       // success
                        }else{  // s2!=sl && s2!=ZERO
                                CGALRS_dyadic_clear(xkappa);
                                CGALRS_dyadic_clear(w);
                                CGALRS_dyadic_clear(xside);
                                return -1;      // failure
                        }
                }
        }
        CGAL_assertion_msg(false,"never reached");
        return 2;
}

// applies qir a given number of times
template <class _Gcd_policy>
int qir_n(const Algebraic_1 &a,unsigned t=1){
        typedef _Gcd_policy     Gcd;
        unsigned count=0;
        unsigned n=2;   // N=2^n
        CGALRS_dyadic_t x_lo,x_hi,f_lo,f_hi; // endpoints and their evaluations
        if(a.is_point())
                return 0;
        CGALRS_dyadic_init_set_fr(x_lo,a.left());
        CGALRS_dyadic_init_set_fr(x_hi,a.right());
        CGALRS_dyadic_init(f_lo);
        CGALRS_dyadic_init(f_hi);
        sfpart_1<Gcd>()(a.pol()).eval_dyadic(f_lo,x_lo);
        sfpart_1<Gcd>()(a.pol()).eval_dyadic(f_hi,x_hi);
        if(!CGALRS_dyadic_sgn(f_lo)){
                CGALRS_dyadic_get_exactfr(&(a.mpfi()->right),x_lo);
                return 0;
        }
        if(!CGALRS_dyadic_sgn(f_hi)){
                CGALRS_dyadic_get_exactfr(&(a.mpfi()->left),x_hi);
                return 0;
        }
        CGAL_assertion(CGALRS_dyadic_sgn(f_lo)!=CGALRS_dyadic_sgn(f_hi));
        while(t>count){
                switch(refine_interval_by_factor(x_lo,x_hi,f_lo,f_hi,
                                        sfpart_1<Gcd>()(a.pol()),n))
                {
                        case 0:n*=2;++count;break;      // N=N^2
                        case -1:if(n>2)n/=2;break;      // N=sqrt(N)
                        case -2:++count;break;
                        default:++count;goto endloop;   // we found the root
                }
        }
        endloop:
        CGALRS_dyadic_get_exactfr(&(a.mpfi()->left),x_lo);
        CGALRS_dyadic_get_exactfr(&(a.mpfi()->right),x_hi);
        CGALRS_dyadic_clear(x_lo);
        CGALRS_dyadic_clear(x_hi);
        CGALRS_dyadic_clear(f_lo);
        CGALRS_dyadic_clear(f_hi);

        a.set_lefteval(RSSign::signat(sfpart_1<Gcd>()(a.pol()),a.left()));

        return count;
}

// applies qir until having an interval of size less than 2^(-t)
template <class _Gcd_policy>
int qir(const Algebraic_1 &a,mp_exp_t t){
        typedef _Gcd_policy     Gcd;
        int count=0;
        long ed;
        unsigned n=2;   // N=2^n
        CGALRS_dyadic_t d,f_lo,f_hi;
        if(a.is_point())
                return 0;
        CGALRS_dyadic_init(f_lo);
        CGALRS_dyadic_init(f_hi);
        sfpart_1<Gcd>()
                (a.pol()).eval_dyadic(f_lo,(CGALRS_dyadic_ptr)(a.left()));
        sfpart_1<Gcd>()
                (a.pol()).eval_dyadic(f_hi,(CGALRS_dyadic_ptr)(a.right()));
        // we make sure we have a nice interval
        if(!CGALRS_dyadic_sgn(f_lo)){
                CGALRS_dyadic_set((CGALRS_dyadic_ptr)a.right(),
                                  (CGALRS_dyadic_srcptr)a.left());
                return 0;
        }
        if(!CGALRS_dyadic_sgn(f_hi)){
                CGALRS_dyadic_set((CGALRS_dyadic_ptr)a.left(),
                                  (CGALRS_dyadic_srcptr)a.right());
                return 0;
        }
        CGAL_assertion(CGALRS_dyadic_sgn(f_lo)!=CGALRS_dyadic_sgn(f_hi));
        // calculate ed, such that 2^(ed-1)<diam<2^ed
        CGALRS_dyadic_init(d);
        CGALRS_dyadic_sub(d,
                          (CGALRS_dyadic_ptr)a.right(),
                          (CGALRS_dyadic_ptr)a.left());
        mpfr_log2((mpfr_ptr)d,(mpfr_ptr)d,GMP_RNDU);
        ed=mpfr_get_si((mpfr_ptr)d,GMP_RNDU);
        while(ed<t){
                ++count;
                switch(refine_interval_by_factor(
                                (CGALRS_dyadic_ptr)a.left(),
                                (CGALRS_dyadic_ptr)a.right(),
                                f_lo,
                                f_hi,
                                sfpart_1<Gcd>()(a.pol()),
                                n)){
                        case 0:t-=n;n*=2;break; // N=N^2
                        case -1:if(n>2)n/=2;break;      // N=sqrt(N)
                        case -2:t-=2;break;
                        default:t=0;    // end the loop, we found the root
                }
        }
        CGAL_assertion(CGALRS_dyadic_sgn(f_lo)!=CGALRS_dyadic_sgn(f_hi));
        a.set_lefteval(
                CGALRS_dyadic_sgn(f_lo)==0?
                ZERO:
                (CGALRS_dyadic_sgn(f_lo)<0?NEGATIVE:POSITIVE));
        CGALRS_dyadic_clear(d);
        CGALRS_dyadic_clear(f_lo);
        CGALRS_dyadic_clear(f_hi);

        return count;
}

} // namespace CGAL

#endif  // CGAL_RS_REFINE_1_H

// vim: tabstop=8: softtabstop=8: smarttab: shiftwidth=8: expandtab
