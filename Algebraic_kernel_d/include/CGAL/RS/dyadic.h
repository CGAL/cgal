// Copyright (c) 2006-2013 INRIA Nancy-Grand Est (France). All rights reserved.
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
// Author: Luis Peñaranda <luis.penaranda@gmx.com>

#ifndef CGAL_RS_DYADIC_H
#define CGAL_RS_DYADIC_H

#include <stdio.h>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include <CGAL/assertions.h>

// for c++, compile with -lgmpxx
#ifdef __cplusplus
#include <iostream>
#endif

#define CGALRS_dyadic_struct            __mpfr_struct
#define CGALRS_dyadic_t                 mpfr_t
#define CGALRS_dyadic_ptr               mpfr_ptr
#define CGALRS_dyadic_srcptr            mpfr_srcptr

// some auxiliary defines
#define CGALRS_dyadic_set_prec(D,P) \
 ( mpfr_set_prec( (D), (P)>MPFR_PREC_MIN?(P):MPFR_PREC_MIN) )

#define CGALRS_dyadic_prec_round(D,P) \
 ( mpfr_prec_round( (D), (P)>MPFR_PREC_MIN?(P):MPFR_PREC_MIN, GMP_RNDN) )

#define CGALRS_dyadic_set_exp(D,E) \
 ( CGAL_assertion( (E) <= mpfr_get_emax() && \
                   (E) >= mpfr_get_emin() ) ,\
   mpfr_set_exp(D,E) )

// init functions
#define CGALRS_dyadic_init(D)          mpfr_init2(D,MPFR_PREC_MIN)
#define CGALRS_dyadic_init2(D,P)       mpfr_init2(D,P)
#define CGALRS_dyadic_clear(D)         mpfr_clear(D)

inline void CGALRS_dyadic_set(CGALRS_dyadic_ptr rop,CGALRS_dyadic_srcptr op){
        if(rop!=op){
                CGALRS_dyadic_set_prec(rop,mpfr_get_prec(op));
                mpfr_set(rop,op,GMP_RNDN);
        }
        CGAL_assertion(mpfr_equal_p(rop,op)!=0);
}

inline void CGALRS_dyadic_set_z(CGALRS_dyadic_ptr rop,mpz_srcptr z){
        size_t prec;
        prec=mpz_sizeinbase(z,2)-(mpz_tstbit(z,0)?0:mpz_scan1(z,0));
        CGALRS_dyadic_set_prec(rop,prec);
        mpfr_set_z(rop,z,GMP_RNDN);
        CGAL_assertion(!mpfr_cmp_z(rop,z));
}

inline void CGALRS_dyadic_set_si(CGALRS_dyadic_ptr rop,long s){
        CGALRS_dyadic_set_prec(rop,sizeof(long));
        mpfr_set_si(rop,s,GMP_RNDN);
        CGAL_assertion(!mpfr_cmp_si(rop,s));
}

inline void CGALRS_dyadic_set_ui(CGALRS_dyadic_ptr rop,unsigned long u){
        CGALRS_dyadic_set_prec(rop,sizeof(unsigned long));
        mpfr_set_ui(rop,u,GMP_RNDN);
        CGAL_assertion(!mpfr_cmp_ui(rop,u));
}

inline void CGALRS_dyadic_set_fr(CGALRS_dyadic_ptr rop,mpfr_srcptr op){
        if(rop!=op){
                CGALRS_dyadic_set_prec(rop,mpfr_get_prec(op));
                mpfr_set(rop,op,GMP_RNDN);
                CGAL_assertion(mpfr_equal_p(rop,op)!=0);
        }
}

#define CGALRS_dyadic_init_set(R,D) \
 ( CGALRS_dyadic_init(R), CGALRS_dyadic_set((R), (D)) )
#define CGALRS_dyadic_init_set_z(R,Z) \
 ( CGALRS_dyadic_init(R), CGALRS_dyadic_set_z((R), (Z)) )
#define CGALRS_dyadic_init_set_si(R,I) \
 ( CGALRS_dyadic_init(R), CGALRS_dyadic_set_si((R), (I)) )
#define CGALRS_dyadic_init_set_ui(R,I) \
 ( CGALRS_dyadic_init(R), CGALRS_dyadic_set_ui((R), (I)) )
#define CGALRS_dyadic_init_set_fr(R,F) \
 ( CGALRS_dyadic_init(R), CGALRS_dyadic_set_fr((R), (F)) )

#define CGALRS_dyadic_get_fr(M,D)      mpfr_set(M,D,GMP_RNDN)
#define CGALRS_dyadic_get_d(D,RM)      mpfr_get_d(D,RM)
inline void CGALRS_dyadic_get_exactfr(mpfr_ptr rop,CGALRS_dyadic_srcptr op){
        if(rop!=op){
                CGALRS_dyadic_set_prec(rop,mpfr_get_prec(op));
                mpfr_set(rop,op,GMP_RNDN);
                CGAL_assertion(mpfr_equal_p(rop,op)!=0);
        }
}

#define CGALRS_dyadic_canonicalize(D)  ()

// comparison functions
#define CGALRS_dyadic_sgn(D)    mpfr_sgn(D)
#define CGALRS_dyadic_zero(D)   mpfr_zero_p(D)
#define CGALRS_dyadic_cmp(D,E)  mpfr_cmp(D,E)

// arithmetic functions
#define CGALRS_dyadic_add(R,D,E)        CGALRS_dyadic_ll_add(R,D,E,0)
#define CGALRS_dyadic_sub(R,D,E)        CGALRS_dyadic_ll_sub(R,D,E,0)
#define CGALRS_dyadic_mul(R,D,E)        CGALRS_dyadic_ll_mul(R,D,E,0)

inline void CGALRS_dyadic_neg(CGALRS_dyadic_ptr rop,CGALRS_dyadic_srcptr op){
        if(rop!=op)
                CGALRS_dyadic_set_prec(rop,mpfr_get_prec(op));
        mpfr_neg(rop,op,GMP_RNDN);
        CGAL_assertion(
                rop==op||
                (!mpfr_cmpabs(rop,op)&&
                ((CGALRS_dyadic_zero(op)&&CGALRS_dyadic_zero(rop))||
                 (CGALRS_dyadic_sgn(op)!=CGALRS_dyadic_sgn(rop)))));
}

// low-level addition:
// add op1 and op2 and reserve b bits for future lowlevel operations
inline void CGALRS_dyadic_ll_add(CGALRS_dyadic_ptr rop,
                                 CGALRS_dyadic_srcptr op1,
                                 CGALRS_dyadic_srcptr op2,
                                 mp_prec_t b){
        mp_exp_t l,r,temp1,temp2;
        mp_prec_t rop_prec;
        if(mpfr_zero_p(op1)){
                if(rop!=op2)
                        CGALRS_dyadic_set(rop,op2);
                return;
        }
        if(mpfr_zero_p(op2)){
                if(rop!=op1)
                        CGALRS_dyadic_set(rop,op1);
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
                CGALRS_dyadic_prec_round(rop,rop_prec);
        else
                CGALRS_dyadic_set_prec(rop,rop_prec);
        CGAL_assertion_code(int round=)
        mpfr_add(rop,op1,op2,GMP_RNDN);
        CGAL_assertion(!round);
}

inline void CGALRS_dyadic_add_z(CGALRS_dyadic_ptr rop,
                                CGALRS_dyadic_srcptr op1,
                                mpz_srcptr z){
        mp_exp_t l,r;
        mp_prec_t rop_prec;
        if(mpfr_zero_p(op1)){
                CGALRS_dyadic_set_z(rop,z);
                return;
        }
        if(!mpz_sgn(z)){
                if(rop!=op1)
                        CGALRS_dyadic_set(rop,op1);
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
                CGALRS_dyadic_prec_round(rop,rop_prec);
        else
                CGALRS_dyadic_set_prec(rop,rop_prec);
        CGAL_assertion_code(int round=)
        mpfr_add_z(rop,op1,z,GMP_RNDN);
        CGAL_assertion(!round);
}

// low-level subtraction:
// subtract op2 to op1 and reserve b bits for future lowlevel operations
inline void CGALRS_dyadic_ll_sub(CGALRS_dyadic_ptr rop,
                                 CGALRS_dyadic_srcptr op1,
                                 CGALRS_dyadic_srcptr op2,
                                 mp_prec_t b){
        mp_exp_t l,r,temp1,temp2;
        mp_prec_t rop_prec;
        if(mpfr_zero_p(op1)){
                CGALRS_dyadic_neg(rop,op2);
                return;
        }
        if(mpfr_zero_p(op2)){
                if(rop!=op1)
                        CGALRS_dyadic_set(rop,op1);
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
                CGALRS_dyadic_prec_round(rop,rop_prec);
        else
                CGALRS_dyadic_set_prec(rop,rop_prec);
        CGAL_assertion_code(int round=)
        mpfr_sub(rop,op1,op2,GMP_RNDN);
        CGAL_assertion(!round);
}

// low-level multiplication:
// multiply op1 and op2 and reserve b bits for future lowlevel operations
inline void CGALRS_dyadic_ll_mul(CGALRS_dyadic_ptr rop,
                                 CGALRS_dyadic_srcptr op1,
                                 CGALRS_dyadic_srcptr op2,
                                 mp_prec_t b){
        if(rop==op1||rop==op2)
                CGALRS_dyadic_prec_round(rop,b+mpfr_get_prec(op1)+mpfr_get_prec(op2));
        else
                CGALRS_dyadic_set_prec(rop,b+mpfr_get_prec(op1)+mpfr_get_prec(op2));
        CGAL_assertion_code(int round=)
        mpfr_mul(rop,op1,op2,GMP_RNDN);
        CGAL_assertion(!round);
}

inline void CGALRS_dyadic_mul_z(CGALRS_dyadic_ptr rop,
                                CGALRS_dyadic_srcptr op1,
                                mpz_srcptr z){
        if(rop==op1)
                CGALRS_dyadic_prec_round(
                        rop,
                        mpfr_get_prec(op1)+mpz_sizeinbase(z,2));
        else
                CGALRS_dyadic_set_prec(
                        rop,
                        mpfr_get_prec(op1)+mpz_sizeinbase(z,2));
        CGAL_assertion_code(int round=)
        mpfr_mul_z(rop,op1,z,GMP_RNDN);
        CGAL_assertion(!round);
}

inline void CGALRS_dyadic_mul_si(CGALRS_dyadic_ptr rop,
                                 CGALRS_dyadic_srcptr op1,
                                 long s){
        if(rop==op1)
                CGALRS_dyadic_prec_round(rop,mpfr_get_prec(op1)+sizeof(long));
        else
                CGALRS_dyadic_set_prec(rop,mpfr_get_prec(op1)+sizeof(long));
        CGAL_assertion_code(int round=)
        mpfr_mul_si(rop,op1,s,GMP_RNDN);
        CGAL_assertion(!round);
}

inline void CGALRS_dyadic_mul_ui(CGALRS_dyadic_ptr rop,
                                 CGALRS_dyadic_srcptr op1,
                                 unsigned long u){
        if(rop==op1)
                CGALRS_dyadic_prec_round(
                        rop,
                        mpfr_get_prec(op1)+sizeof(unsigned long));
        else
                CGALRS_dyadic_set_prec(
                        rop,
                        mpfr_get_prec(op1)+sizeof(unsigned long));
        CGAL_assertion_code(int round=)
        mpfr_mul_ui(rop,op1,u,GMP_RNDN);
        CGAL_assertion(!round);
}

inline void CGALRS_dyadic_pow_ui(CGALRS_dyadic_ptr rop,
                                 CGALRS_dyadic_srcptr op1,
                                 unsigned long u){
        if(!u){
                CGAL_assertion_msg(!mpfr_zero_p(op1),"0^0");
                CGALRS_dyadic_set_ui(rop,1);
                return;
        }
        if(u==1){
                if(rop!=op1)
                        CGALRS_dyadic_set(rop,op1);
                return;
        }
        if(mpfr_zero_p(op1)){
                CGAL_assertion_msg(u!=0,"0^0");
                CGALRS_dyadic_set_ui(rop,0);
                return;
        }
        if(!mpfr_cmp_ui(op1,1)){
                if(rop!=op1)
                        CGALRS_dyadic_set(rop,op1);
                return;
        }
        if(rop==op1)
                CGALRS_dyadic_prec_round(rop,u*mpfr_get_prec(op1));
        else
                CGALRS_dyadic_set_prec(rop,u*mpfr_get_prec(op1));
        CGAL_assertion_code(int round=)
        mpfr_pow_ui(rop,op1,u,GMP_RNDN);
        CGAL_assertion(!round);
}

inline void CGALRS_dyadic_addmul(CGALRS_dyadic_ptr rop,
                                 CGALRS_dyadic_srcptr op1,
                                 CGALRS_dyadic_srcptr op2){
        CGALRS_dyadic_t temp;
        CGALRS_dyadic_init(temp);
        CGALRS_dyadic_mul(temp,op1,op2);
        CGALRS_dyadic_add(rop,rop,temp);
        CGALRS_dyadic_clear(temp);
}

inline void CGALRS_dyadic_addmul_si(CGALRS_dyadic_ptr rop,
                                    CGALRS_dyadic_srcptr op1,
                                    long op2){
        CGALRS_dyadic_t temp;
        CGALRS_dyadic_init(temp);
        CGALRS_dyadic_mul_si(temp,op1,op2);
        CGALRS_dyadic_add(rop,rop,temp);
        CGALRS_dyadic_clear(temp);
}

inline void CGALRS_dyadic_addmul_ui(CGALRS_dyadic_ptr rop,
                                    CGALRS_dyadic_srcptr op1,
                                    unsigned long u){
        //CGALRS_dyadic_t temp;
        //CGALRS_dyadic_init(temp);
        //CGALRS_dyadic_mul_ui(temp,op1,u);
        //CGALRS_dyadic_add(rop,rop,temp);
        //CGALRS_dyadic_clear(temp);
        CGALRS_dyadic_t temp;
        mp_exp_t l,r,temp1,temp2;
        mp_prec_t rop_prec;
        if(u==0||mpfr_zero_p(op1))
                return;
        if(u==1){
                CGALRS_dyadic_add(rop,rop,op1);
                return;
        }
        // TODO: if(op1==1)
        // calculate temp=op1*u
        mpfr_init2(temp,mpfr_get_prec(op1)+sizeof(unsigned int));
        CGAL_assertion_code(int round1=)
        mpfr_mul_ui(temp,op1,u,GMP_RNDN);
        CGAL_assertion(!round1);
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
        CGALRS_dyadic_prec_round(rop,rop_prec);
        CGAL_assertion_code(int round2=)
        mpfr_add(rop,rop,temp,GMP_RNDN);
        CGAL_assertion(!round2);
}

inline void CGALRS_dyadic_mul_2exp(CGALRS_dyadic_ptr rop,
                                   CGALRS_dyadic_srcptr op1,
                                   unsigned long ui){
        // mpfr_mul_2ui does change the mantissa!
        if(rop==op1)
                CGALRS_dyadic_prec_round(
                        rop,
                        sizeof(unsigned long)+mpfr_get_prec(op1));
        else
                CGALRS_dyadic_set_prec(
                        rop,
                        sizeof(unsigned long)+mpfr_get_prec(op1));
        CGAL_assertion_code(int round=)
        mpfr_mul_2ui(rop,op1,ui,GMP_RNDN);
        CGAL_assertion(!round);
}

inline void CGALRS_dyadic_div_2exp(CGALRS_dyadic_ptr rop,
                                   CGALRS_dyadic_srcptr op1,
                                   unsigned long ui){
        // mpfr_div_2ui does not change the mantissa... am I sure?
        CGAL_assertion_code(int round=)
        mpfr_div_2ui(rop,op1,ui,GMP_RNDN);
        CGAL_assertion(!round);
}

// miscelaneous functions
#define CGALRS_dyadic_midpoint(R,D,E) \
 ( CGALRS_dyadic_ll_add(R,D,E,1) , mpfr_div_2ui(R,R,1,GMP_RNDN) )
#define CGALRS_dyadic_swap(D,E)         mpfr_swap(D,E)

// I/O functions
#define CGALRS_dyadic_out_str(F,D)      mpfr_out_str(F,10,0,D,GMP_RNDN)
#ifdef __cplusplus
inline std::ostream& operator<<(std::ostream &s,CGALRS_dyadic_srcptr op){
        mp_exp_t exponent;
        mpz_t mantissa;
        mpz_init(mantissa);
        exponent=mpfr_get_z_exp(mantissa,op);
        s<<"["<<mantissa<<","<<exponent<<"]";
        return s;
}
#endif

#endif  // CGAL_RS_DYADIC_H
