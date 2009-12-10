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

#ifndef CGAL_RS__P_H
#define CGAL_RS__P_H

#include <CGAL/RS/basic.h>
#include <CGAL/RS/inverse.h>

namespace CGAL{
namespace RS_MGCD{

// pn size is 32 bits,
// the sizes of lpn and spn must be, at least, as twice as the size of pn
#define PN_BITS 32
#define pn      uint32_t        // unsigned
#define lpn     uint64_t        // unsigned long long
#define spn     int64_t         // long long

#define p_mpz_set_pn(A,PN)      mpz_set_ui(A,(unsigned long)(PN))
#define p_mpz_mul_pn(A,B,PN)    mpz_mul_ui(A,B,(unsigned long)(PN))
#define p_mpz_add_pn(A,B,PN)    mpz_add_ui(A,B,(unsigned long)(PN))
#define p_mpz_sub_pn(A,B,PN)    mpz_sub_ui(A,B,(unsigned long)(PN))
#define p_mpz_set_spn(A,SPN)    mpz_set_si(A,(long)(SPN))
#define p_mpz_add_spn(A,B,SPN)  (SPN<0?p_mpz_sub_pn(A,B,-(SPN)):p_mpz_add_pn(A,B,(SPN)))
#define p_mpz_mul_spn(A,B,SPN)  mpz_mul_si(A,B,SPN)

RS_THREAD_ATTR pn prime;

class Prime:public Inverse{
    protected:
        static spn p_pntospn(pn p){
            if(p>(prime-1)/2)
                return (spn)p-prime;
            return (spn)p;
        };

        static void p_set_prime(pn p){prime=p;};

        static pn p_prime(){return prime;};

        static pn p_add(pn a,pn b){
            lpn c=(lpn)a+b;
            return (c<prime?(pn)c:(pn)(c-prime));
        };

        static pn p_sub(pn a,pn b){
            return (a<b?prime-(b-a):a-b);
        };

        // returns a*b-c*d
        static pn p_submuls(pn a,pn b,pn c,pn d){
            lpn mul;
            pn pnm1,pnm2;
            mul=(lpn)a*b;
            pnm1=(mul<prime?(pn)mul:(pn)(mul%prime));
            mul=(lpn)c*d;
            pnm2=(mul<prime?(pn)mul:(pn)(mul%prime));
            return (pnm1<pnm2?prime-(pnm2-pnm1):pnm1-pnm2);
        };

        static pn p_inv(pn a){
            spn c=eea_s(a,prime);
            return(c<0?(pn)(c+prime):(pn)c);
        };

        static pn p_convert(pn a){
            return (a<prime?a:a%prime);
        };

        //#define p_mul(A,B)    ((pn)(((lpn)A*B)%p_prime()))
        static pn p_mul(pn a,pn b){
            lpn c=(lpn)a*b;
            return (pn)(c<prime?c:c%prime);
        };

        static pn p_mul3(pn a,pn b,pn c){
            lpn d;
            if((d=(lpn)a*b)<prime)
                d*=c;
            else
                d=(d%prime)*c;
            return (pn)(d<prime?d:d%prime);
        };

        // returns a*conv(b)
        static pn p_mulc(pn a,pn b){
            lpn temp=(lpn)a*(b<prime?b:b%prime);
            return (temp<prime?(pn)temp:(pn)(temp%prime));
        };

        // returns a*conv(b)+conv(c)
        static pn p_mulcaddc(pn a,pn b,pn c){
            lpn temp=(lpn)a*(b<prime?b:b%prime);
            temp=(temp<prime?temp:temp%prime)+(c<prime?c:c%prime);
            return (temp<prime?(pn)temp:(pn)(temp-prime));
        };

        // returns (conv(a)-b)*inv(c)
        static pn p_convsubdiv(pn a,pn b,pn c){
            spn inv_c=eea_s(c,prime);
            pn pninv_c=(inv_c<0?(pn)(inv_c+prime):(pn)inv_c);
            pn conv_a=(a<prime?a:a%prime);
            lpn mult=(lpn)(conv_a<b?prime-(b-conv_a):conv_a-b)*pninv_c;
            return (pn)(mult<prime?mult:mult%prime);
        };

        // vzGG, p. 73
        static pn p_pow(pn a,pn n){
            pn b,i;
            i=1<<(PN_BITS-1);
            while(!(i&n))
                i>>=1;
            b=a;
            while(i>>=1)
                b=(i&n?p_mul3(b,b,a):p_mul(b,b));
            return b;
        };

        #define p_div(A,B)      (p_mul(A,p_inv(B)))

        static pn p_gcd(pn a,pn b){
            if(!b)
                return a;
            return p_gcd(b,a%b);
        };

}; // class Prime

} // namespace RS_MGCD
} // namespace CGAL

#endif  // CGAL_RS__P_H

// vim: tabstop=4: softtabstop=4: smarttab: shiftwidth=4: expandtab
