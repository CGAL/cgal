// Copyright (c) 2007-2008 Inria Lorraine (France). All rights reserved.
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

#ifndef _dyadic_h
#define _dyadic_h

#include <stdio.h>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include <CGAL/assertions.h>

// for c++, compile with -lgmpxx
#ifdef __cplusplus
#include <gmpxx.h>
#include <iostream>
#endif

#define __dyadic_struct __mpfr_struct
#define dyadic_t        mpfr_t
#define dyadic_ptr      mpfr_ptr
#define dyadic_srcptr   mpfr_srcptr

// some auxiliary defines
#define dyadic_set_prec(D,P) \
 ( mpfr_set_prec( (D), (P)>MPFR_PREC_MIN?(P):MPFR_PREC_MIN) )

#define dyadic_prec_round(D,P) \
 ( mpfr_prec_round( (D), (P)>MPFR_PREC_MIN?(P):MPFR_PREC_MIN, GMP_RNDN) )

#define dyadic_set_exp(D,E) \
 ( CGAL_assertion( (E) <= mpfr_get_emax() && \
                   (E) >= mpfr_get_emin() ) ,\
   mpfr_set_exp(D,E) )

// init functions
#define dyadic_init(D)          mpfr_init2(D,MPFR_PREC_MIN)
#define dyadic_init2(D,P)       mpfr_init2(D,P)
#define dyadic_clear(D)         mpfr_clear(D)

inline void dyadic_set(dyadic_ptr rop,dyadic_srcptr op){
        if(rop!=op){
                dyadic_set_prec(rop,mpfr_get_prec(op));
                mpfr_set(rop,op,GMP_RNDN);
        }
        CGAL_assertion(mpfr_equal_p(rop,op));
}

inline void dyadic_set_z(dyadic_ptr rop,mpz_srcptr z){
        size_t prec;
        prec=mpz_sizeinbase(z,2)-(mpz_tstbit(z,0)?0:mpz_scan1(z,0));
        dyadic_set_prec(rop,prec);
        mpfr_set_z(rop,z,GMP_RNDN);
        CGAL_assertion(!mpfr_cmp_z(rop,z));
}

inline void dyadic_set_si(dyadic_ptr rop,long s){
        dyadic_set_prec(rop,sizeof(long));
        mpfr_set_si(rop,s,GMP_RNDN);
        CGAL_assertion(!mpfr_cmp_si(rop,s));
}

inline void dyadic_set_ui(dyadic_ptr rop,unsigned long u){
        dyadic_set_prec(rop,sizeof(unsigned long));
        mpfr_set_ui(rop,u,GMP_RNDN);
        CGAL_assertion(!mpfr_cmp_ui(rop,u));
}

inline void dyadic_set_fr(dyadic_ptr rop,mpfr_srcptr op){
        if(rop!=op){
                dyadic_set_prec(rop,mpfr_get_prec(op));
                mpfr_set(rop,op,GMP_RNDN);
                CGAL_assertion(mpfr_equal_p(rop,op));
        }
}

#define dyadic_init_set(R,D) \
 ( dyadic_init(R), dyadic_set((R), (D)) )
#define dyadic_init_set_z(R,Z) \
 ( dyadic_init(R), dyadic_set_z((R), (Z)) )
#define dyadic_init_set_si(R,I) \
 ( dyadic_init(R), dyadic_set_si((R), (I)) )
#define dyadic_init_set_ui(R,I) \
 ( dyadic_init(R), dyadic_set_ui((R), (I)) )
#define dyadic_init_set_fr(R,F) \
 ( dyadic_init(R), dyadic_set_fr((R), (F)) )

#define dyadic_get_fr(M,D)      mpfr_set(M,D,GMP_RNDN)
#define dyadic_get_d(D,RM)      mpfr_get_d(D,RM)
inline void dyadic_get_exactfr(mpfr_ptr rop,dyadic_srcptr op){
        if(rop!=op){
                dyadic_set_prec(rop,mpfr_get_prec(op));
                mpfr_set(rop,op,GMP_RNDN);
                CGAL_assertion(mpfr_equal_p(rop,op));
        }
}

#define dyadic_canonicalize(D)  ()

// comparison functions
#define dyadic_sgn(D)   mpfr_sgn(D)
#define dyadic_zero(D)  mpfr_zero_p(D)
#define dyadic_cmp(D,E) mpfr_cmp(D,E)

// arithmetic functions
#define dyadic_add(R,D,E)       dyadic_ll_add(R,D,E,0)
#define dyadic_sub(R,D,E)       dyadic_ll_sub(R,D,E,0)
#define dyadic_mul(R,D,E)       dyadic_ll_mul(R,D,E,0)

inline void dyadic_neg(dyadic_ptr rop,dyadic_srcptr op){
        if(rop!=op)
                dyadic_set_prec(rop,mpfr_get_prec(op));
        mpfr_neg(rop,op,GMP_RNDN);
        CGAL_assertion(
                rop==op||
                (!mpfr_cmpabs(rop,op)&&
                ((dyadic_zero(op)&&dyadic_zero(rop))||
                 (dyadic_sgn(op)!=dyadic_sgn(rop)))));
}

// low-level addition:
// add op1 and op2 and reserve b bits for future lowlevel operations
inline void dyadic_ll_add(
                        dyadic_ptr rop,
                        dyadic_srcptr op1,
                        dyadic_srcptr op2,
                        mp_prec_t b){
        mp_exp_t l,r,temp1,temp2;
        mp_prec_t rop_prec;
        int round;
        if(mpfr_zero_p(op1)){
                if(rop!=op2)
                        dyadic_set(rop,op2);
                return;
        }
        if(mpfr_zero_p(op2)){
                if(rop!=op1)
                        dyadic_set(rop,op1);
                return;
        }
        l=mpfr_get_exp(op1)>mpfr_get_exp(op2)?
                mpfr_get_exp(op1):
                mpfr_get_exp(op2);
        temp1=mpfr_get_exp(op1)-(mp_exp_t)mpfr_get_prec(op1);
        temp2=mpfr_get_exp(op2)-(mp_exp_t)mpfr_get_prec(op2);
        r=temp1>temp2?temp2:temp1;
        CGAL_assertion(l>r);

        rop_prec=b+1+(mp_prec_t)(l-r);
        CGAL_assertion(rop_prec>=mpfr_get_prec(op1)&&
                rop_prec>=mpfr_get_prec(op2));
        if(rop==op1||rop==op2)
                dyadic_prec_round(rop,rop_prec);
        else
                dyadic_set_prec(rop,rop_prec);
        round=mpfr_add(rop,op1,op2,GMP_RNDN);
        CGAL_assertion(!round);
}

inline void dyadic_add_z(dyadic_ptr rop,dyadic_srcptr op1,mpz_srcptr z){
        mp_exp_t l,r;
        mp_prec_t rop_prec;
        int round;
        if(mpfr_zero_p(op1)){
                dyadic_set_z(rop,z);
                return;
        }
        if(!mpz_sgn(z)){
                if(rop!=op1)
                        dyadic_set(rop,op1);
                return;
        }
        l=mpfr_get_exp(op1)>(mp_exp_t)mpz_sizeinbase(z,2)?
                mpfr_get_exp(op1):
                mpz_sizeinbase(z,2);
        r=mpfr_get_exp(op1)-(mp_exp_t)mpfr_get_prec(op1)<0?
                mpfr_get_exp(op1)-(mp_exp_t)mpfr_get_prec(op1):
                0;
        CGAL_assertion(l>r);

        rop_prec=1+(mp_prec_t)(l-r);
        CGAL_assertion(rop_prec>=mpfr_get_prec(op1)&&
                rop_prec>=(mp_prec_t)mpz_sizeinbase(z,2));
        if(rop==op1)
                dyadic_prec_round(rop,rop_prec);
        else
                dyadic_set_prec(rop,rop_prec);
        round=mpfr_add_z(rop,op1,z,GMP_RNDN);
        CGAL_assertion(!round);
}

// low-level subtraction:
// subtract op2 to op1 and reserve b bits for future lowlevel operations
inline void dyadic_ll_sub(
                        dyadic_ptr rop,
                        dyadic_srcptr op1,
                        dyadic_srcptr op2,
                        mp_prec_t b){
        mp_exp_t l,r,temp1,temp2;
        mp_prec_t rop_prec;
        int round;
        if(mpfr_zero_p(op1)){
                dyadic_neg(rop,op2);
                return;
        }
        if(mpfr_zero_p(op2)){
                if(rop!=op1)
                        dyadic_set(rop,op1);
                return;
        }
        l=mpfr_get_exp(op1)>mpfr_get_exp(op2)?
                mpfr_get_exp(op1):
                mpfr_get_exp(op2);
        temp1=mpfr_get_exp(op1)-(mp_exp_t)mpfr_get_prec(op1);
        temp2=mpfr_get_exp(op2)-(mp_exp_t)mpfr_get_prec(op2);
        r=temp1>temp2?temp2:temp1;
        CGAL_assertion(l>r);

        rop_prec=b+1+(mp_prec_t)(l-r);
        CGAL_assertion(rop_prec>=mpfr_get_prec(op1)&&
                rop_prec>=mpfr_get_prec(op2));
        if(rop==op1||rop==op2)
                dyadic_prec_round(rop,rop_prec);
        else
                dyadic_set_prec(rop,rop_prec);
        round=mpfr_sub(rop,op1,op2,GMP_RNDN);
        CGAL_assertion(!round);
}

// low-level multiplication:
// multiply op1 and op2 and reserve b bits for future lowlevel operations
inline void dyadic_ll_mul(
                        dyadic_ptr rop,
                        dyadic_srcptr op1,
                        dyadic_srcptr op2,
                        mp_prec_t b){
        int round;
        if(rop==op1||rop==op2)
                dyadic_prec_round(rop,b+mpfr_get_prec(op1)+mpfr_get_prec(op2));
        else
                dyadic_set_prec(rop,b+mpfr_get_prec(op1)+mpfr_get_prec(op2));
        round=mpfr_mul(rop,op1,op2,GMP_RNDN);
        CGAL_assertion(!round);
}

inline void dyadic_mul_z(dyadic_ptr rop,dyadic_srcptr op1,mpz_srcptr z){
        int round;
        if(rop==op1)
                dyadic_prec_round(rop,mpfr_get_prec(op1)+mpz_sizeinbase(z,2));
        else
                dyadic_set_prec(rop,mpfr_get_prec(op1)+mpz_sizeinbase(z,2));
        round=mpfr_mul_z(rop,op1,z,GMP_RNDN);
        CGAL_assertion(!round);
}

inline void dyadic_mul_si(dyadic_ptr rop,dyadic_srcptr op1,long s){
        int round;
        if(rop==op1)
                dyadic_prec_round(rop,mpfr_get_prec(op1)+sizeof(long));
        else
                dyadic_set_prec(rop,mpfr_get_prec(op1)+sizeof(long));
        round=mpfr_mul_si(rop,op1,s,GMP_RNDN);
        CGAL_assertion(!round);
}

inline void dyadic_mul_ui(dyadic_ptr rop,dyadic_srcptr op1,unsigned long u){
        int round;
        if(rop==op1)
                dyadic_prec_round(rop,mpfr_get_prec(op1)+sizeof(unsigned long));
        else
                dyadic_set_prec(rop,mpfr_get_prec(op1)+sizeof(unsigned long));
        round=mpfr_mul_ui(rop,op1,u,GMP_RNDN);
        CGAL_assertion(!round);
}

inline void dyadic_pow_ui(dyadic_ptr rop,dyadic_srcptr op1,unsigned long u){
        int round;
        if(!u){
                CGAL_assertion_msg(!mpfr_zero_p(op1),"0^0");
                dyadic_set_ui(rop,1);
                return;
        }
        if(u==1){
                if(rop!=op1)
                        dyadic_set(rop,op1);
                return;
        }
        if(mpfr_zero_p(op1)){
                CGAL_assertion_msg(u,"0^0");
                dyadic_set_ui(rop,0);
                return;
        }
        if(!mpfr_cmp_ui(op1,1)){
                if(rop!=op1)
                        dyadic_set(rop,op1);
                return;
        }
        if(rop==op1)
                dyadic_prec_round(rop,u*mpfr_get_prec(op1));
        else
                dyadic_set_prec(rop,u*mpfr_get_prec(op1));
        round=mpfr_pow_ui(rop,op1,u,GMP_RNDN);
        CGAL_assertion(!round);
}

inline void dyadic_addmul(dyadic_ptr rop,dyadic_srcptr op1,dyadic_srcptr op2){
        dyadic_t temp;
        dyadic_init(temp);
        dyadic_mul(temp,op1,op2);
        dyadic_add(rop,rop,temp);
        dyadic_clear(temp);
}

inline void dyadic_addmul_si(dyadic_ptr rop,dyadic_srcptr op1,long op2){
        dyadic_t temp;
        dyadic_init(temp);
        dyadic_mul_si(temp,op1,op2);
        dyadic_add(rop,rop,temp);
        dyadic_clear(temp);
}

inline void dyadic_addmul_ui(dyadic_ptr rop,dyadic_srcptr op1,unsigned long u){
        //dyadic_t temp;
        //dyadic_init(temp);
        //dyadic_mul_ui(temp,op1,u);
        //dyadic_add(rop,rop,temp);
        //dyadic_clear(temp);
        dyadic_t temp;
        mp_exp_t l,r,temp1,temp2;
        mp_prec_t rop_prec;
        int round;
        if(u==0||mpfr_zero_p(op1))
                return;
        if(u==1){
                dyadic_add(rop,rop,op1);
                return;
        }
        // TODO: if(op1==1)
        // calculate temp=op1*u
        mpfr_init2(temp,mpfr_get_prec(op1)+sizeof(unsigned int));
        round=mpfr_mul_ui(temp,op1,u,GMP_RNDN);
        CGAL_assertion(!round);
        // calculate the precision needed for rop
        l=mpfr_get_exp(op1)>0?mpfr_get_exp(op1):0;
        temp1=mpfr_get_exp(op1)-(mp_exp_t)mpfr_get_prec(op1);
        temp2=sizeof(unsigned long);
        r=temp1>temp2?temp2:temp1;
        CGAL_assertion(l>r);
        rop_prec=sizeof(unsigned long)+1+(mp_prec_t)(l-r);
        CGAL_assertion(rop_prec>=mpfr_get_prec(op1)&&
                rop_prec>=mpfr_get_prec(rop));
        // set precision and add
        dyadic_prec_round(rop,rop_prec);
        round=mpfr_add(rop,rop,temp,GMP_RNDN);
        CGAL_assertion(!round);
}

inline void dyadic_mul_2exp(dyadic_ptr rop,dyadic_srcptr op1,unsigned long ui){
        int round;
        // mpfr_mul_2ui does change the mantissa!
        if(rop==op1)
                dyadic_prec_round(rop,sizeof(unsigned long)+mpfr_get_prec(op1));
        else
                dyadic_set_prec(rop,sizeof(unsigned long)+mpfr_get_prec(op1));
        round=mpfr_mul_2ui(rop,op1,ui,GMP_RNDN);
        CGAL_assertion(!round);
}

inline void dyadic_div_2exp(dyadic_ptr rop,dyadic_srcptr op1,unsigned long ui){
        int round;
        // mpfr_div_2ui does not change the mantissa... am I sure?
        round=mpfr_div_2ui(rop,op1,ui,GMP_RNDN);
        CGAL_assertion(!round);
}

// miscelaneous functions
#define dyadic_midpoint(R,D,E) \
 ( dyadic_ll_add(R,D,E,1) , mpfr_div_2ui(R,R,1,GMP_RNDN) )
#define dyadic_swap(D,E)        mpfr_swap(D,E)

// I/O functions
#define dyadic_out_str(F,D)     mpfr_out_str(F,10,0,D,GMP_RNDN)
#ifdef __cplusplus
inline std::ostream& operator<<(std::ostream &s,dyadic_srcptr op){
        mp_exp_t exponent;
        mpz_t mantissa;
        mpz_init(mantissa);
        exponent=mpfr_get_z_exp(mantissa,op);
        s<<"["<<mantissa<<","<<exponent<<"]";
        return s;
}
#endif

#endif  // _dyadic_h

// vim: tabstop=8: softtabstop=8: smarttab: shiftwidth=8: expandtab
