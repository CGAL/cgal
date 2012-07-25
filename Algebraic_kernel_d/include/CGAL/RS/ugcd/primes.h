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

#ifndef CGAL_RS__PRIMES_H
#define CGAL_RS__PRIMES_H

#include "crt.h"

// I borrowed these numbers from Fabrice, this leaves us around 250000 primes
#define CGALRS_PR_MIN 2145338339
#define CGALRS_PR_MAX 2147483647

//#define CGALRS_PR_IS_PRIME(N)        (pr_fermat(N)?pr_is_prime_bruteforce(N):0)
#define CGALRS_PR_IS_PRIME(N)  pr_mrj(N)

#ifdef _MSC_VER
#  define CGAL_RS_RANDOM  rand
#else
#  define CGAL_RS_RANDOM  random
#endif

namespace CGAL{
namespace RS_MGCD{

CGALRS_THREAD_ATTR CGALRS_PN currentprime;

class Primes:public Crt{
    private:
        static int pr_is_prime_bruteforce(CGALRS_PN n){
            int i,sqrtn;
            sqrtn=(CGALRS_PN)(sqrt((double)n));
            for(i=3;i<=sqrtn;++i)
                if(!(n%i))
                    return 0;
            return 1;
        }

        // vzGG, p. 507; returns 0 if n is composite
        static int pr_fermat(CGALRS_PN n){
            p_set_prime(n);
            return(p_pow(2+((CGALRS_PN)CGAL_RS_RANDOM())%(n-4),n-1)==1);
        }

        // Solovay-Strassen
        static int pr_ss(CGALRS_PN n){
            CGALRS_PN a,x;
            p_set_prime(n);
            a=1+(CGALRS_PN)CGAL_RS_RANDOM()%(n-2);
            x=CGALRS_P_DIV(a,n);
            return(!x||p_pow(a,n>>1)!=x);
        }

        // Miller-Rabin
        static int pr_mr(CGALRS_PN n){
            CGALRS_PN s,d,a,x,r;
            s=1;
            d=(n-1)>>1;
            while(!(d&1)){
                ++s;
                d=d>>1;
            }
            p_set_prime(n);
            a=2+(CGALRS_PN)CGAL_RS_RANDOM()%(n-4);
            x=p_pow(a,d);
            if(x==1||x==n-1)
                return 1; // pobably prime
            for(r=1;r<s;++r){
                x=p_mul(x,x);
                if(x==1)
                    return 0; // composite
                if(x==n-1)
                    return 1; // probably
            }
            return 0; // composite
        }

        // Miller-Rabin-Jaeschke
        static int pr_mrj(CGALRS_PN n){
            CGALRS_PN s,d,a[3],r;//,x;
            int index;
            a[0]=2;
            a[1]=7;
            a[2]=61;
            s=1;
            d=(n-1)>>1;
            while(!(d&1)){
                ++s;
                d=d>>1;
            }
            p_set_prime(n);
            index=-1;
        start_test:
            ++index;
            //if(index=3)
            if(index==3)
                return 1; // prime
            if(p_pow(a[index],d)==1)
                goto start_test;
            for(r=0;r<s;++r)
                if(p_pow(a[index],(d<<r))==n-1)
                    goto start_test;
            return 0; // composite
        }

        static CGALRS_PN pr_actual(){
            return currentprime;
        }

    protected:
        static int pr_init(){
            currentprime=CGALRS_PR_MIN;
            return 0;
        }

        static CGALRS_PN pr_next(){
            do{
                currentprime+=2;
            }while(!CGALRS_PR_IS_PRIME(currentprime));
            return currentprime;
        }

}; // class Primes

} // namespace RS_MGCD
} // namespace CGAL

#endif  // CGAL_RS__PRIMES_H
