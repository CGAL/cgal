// Copyright (c) 2007-2008 Inria Lorraine (France). All rights reserved.
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

#ifndef CGAL_RS__P_H
#define CGAL_RS__P_H

#include <CGAL/RS/basic.h>
#include "inverse.h"

namespace CGAL{
namespace RS_MGCD{

// CGALRS_PN size is 32 bits, the sizes of CGALRS_LPN and CGALRS_SPN must be,
// at least, as twice as the size of CGALRS_PN
#define CGALRS_PN_BITS  32
#define CGALRS_PN       CGALRS_U32  // unsigned
#define CGALRS_LPN      CGALRS_U64  // unsigned long long
#define CGALRS_SPN      CGALRS_S64  // long long

#define CGALRS_mpz_set_pn(A,PN)     mpz_set_ui(A,(unsigned long)(PN))
#define CGALRS_mpz_mul_pn(A,B,PN)   mpz_mul_ui(A,B,(unsigned long)(PN))
#define CGALRS_mpz_add_pn(A,B,PN)   mpz_add_ui(A,B,(unsigned long)(PN))
#define CGALRS_mpz_sub_pn(A,B,PN)   mpz_sub_ui(A,B,(unsigned long)(PN))
#define CGALRS_mpz_set_spn(A,SPN)   mpz_set_si(A,(long)(SPN))
#define CGALRS_mpz_add_spn(A,B,SPN) \
        (SPN<0?CGALRS_mpz_sub_pn(A,B,-(SPN)):CGALRS_mpz_add_pn(A,B,(SPN)))
#define CGALRS_mpz_mul_spn(A,B,SPN) mpz_mul_si(A,B,SPN)

CGALRS_THREAD_ATTR CGALRS_PN prime; class Prime:public Inverse{
    protected:
        static CGALRS_SPN p_pntospn(CGALRS_PN p){
            if(p>(prime-1)/2)
                return (CGALRS_SPN)p-prime;
            return (CGALRS_SPN)p;
        };

        static void p_set_prime(CGALRS_PN p){prime=p;};

        static CGALRS_PN p_prime(){return prime;};

        static CGALRS_PN p_add(CGALRS_PN a,CGALRS_PN b){
            CGALRS_LPN c=(CGALRS_LPN)a+b;
            return (c<prime?(CGALRS_PN)c:(CGALRS_PN)(c-prime));
        };

        static CGALRS_PN p_sub(CGALRS_PN a,CGALRS_PN b){
            return (a<b?prime-(b-a):a-b);
        };

        // returns a*b-c*d
        static CGALRS_PN p_submuls(CGALRS_PN a,CGALRS_PN b,CGALRS_PN c,CGALRS_PN d){
            CGALRS_LPN mul;
            CGALRS_PN pnm1,pnm2;
            mul=(CGALRS_LPN)a*b;
            pnm1=(mul<prime?(CGALRS_PN)mul:(CGALRS_PN)(mul%prime));
            mul=(CGALRS_LPN)c*d;
            pnm2=(mul<prime?(CGALRS_PN)mul:(CGALRS_PN)(mul%prime));
            return (pnm1<pnm2?prime-(pnm2-pnm1):pnm1-pnm2);
        };

        static CGALRS_PN p_inv(CGALRS_PN a){
            CGALRS_SPN c=eea_s(a,prime);
            return(c<0?(CGALRS_PN)(c+prime):(CGALRS_PN)c);
        };

        static CGALRS_PN p_convert(CGALRS_PN a){
            return (a<prime?a:a%prime);
        };

        // define p_mul(A,B)    ((CGALRS_PN)(((CGALRS_LPN)A*B)%p_prime()))
        static CGALRS_PN p_mul(CGALRS_PN a,CGALRS_PN b){
            CGALRS_LPN c=(CGALRS_LPN)a*b;
            return (CGALRS_PN)(c<prime?c:c%prime);
        };

        static CGALRS_PN p_mul3(CGALRS_PN a,CGALRS_PN b,CGALRS_PN c){
            CGALRS_LPN d;
            if((d=(CGALRS_LPN)a*b)<prime)
                d*=c;
            else
                d=(d%prime)*c;
            return (CGALRS_PN)(d<prime?d:d%prime);
        };

        // returns a*conv(b)
        static CGALRS_PN p_mulc(CGALRS_PN a,CGALRS_PN b){
            CGALRS_LPN temp=(CGALRS_LPN)a*(b<prime?b:b%prime);
            return (temp<prime?(CGALRS_PN)temp:(CGALRS_PN)(temp%prime));
        };

        // returns a*conv(b)+conv(c)
        static CGALRS_PN p_mulcaddc(CGALRS_PN a,CGALRS_PN b,CGALRS_PN c){
            CGALRS_LPN temp=(CGALRS_LPN)a*(b<prime?b:b%prime);
            temp=(temp<prime?temp:temp%prime)+(c<prime?c:c%prime);
            return (temp<prime?(CGALRS_PN)temp:(CGALRS_PN)(temp-prime));
        };

        // returns (conv(a)-b)*inv(c)
        static CGALRS_PN p_convsubdiv(CGALRS_PN a,CGALRS_PN b,CGALRS_PN c){
            CGALRS_SPN inv_c=eea_s(c,prime);
            CGALRS_PN pninv_c=(inv_c<0?(CGALRS_PN)(inv_c+prime):(CGALRS_PN)inv_c);
            CGALRS_PN conv_a=(a<prime?a:a%prime);
            CGALRS_LPN mult=
                (CGALRS_LPN)(conv_a<b?prime-(b-conv_a):conv_a-b)*pninv_c;
            return (CGALRS_PN)(mult<prime?mult:mult%prime);
        };

        // vzGG, p. 73
        static CGALRS_PN p_pow(CGALRS_PN a,CGALRS_PN n){
            CGALRS_PN b,i;
            i=1<<(CGALRS_PN_BITS-1);
            while(!(i&n))
                i>>=1;
            b=a;
            while(i>>=1)
                b=(i&n?p_mul3(b,b,a):p_mul(b,b));
            return b;
        };

        #define CGALRS_P_DIV(A,B)      (p_mul(A,p_inv(B)))

        static CGALRS_PN p_gcd(CGALRS_PN a,CGALRS_PN b){
            if(!b)
                return a;
            return p_gcd(b,a%b);
        };

}; // class Prime

} // namespace RS_MGCD
} // namespace CGAL

#endif  // CGAL_RS__P_H
