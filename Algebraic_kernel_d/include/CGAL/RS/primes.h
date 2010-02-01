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

#ifndef CGAL_RS__PRIMES_H
#define CGAL_RS__PRIMES_H

#include <CGAL/RS/crt.h>

// I borrowed these numbers from Fabrice, this leaves us around 250000 primes
#define PR_MIN 2145338339
#define PR_MAX 2147483647

//#define pr_is_prime(N)        (pr_fermat(N)?pr_is_prime_bruteforce(N):0)
#define pr_is_prime(N)  pr_mrj(N)

#ifdef _MSC_VER
#  define cgal_rs_random  rand
#else
#  define cgal_rs_random  random
#endif

namespace CGAL{
namespace RS_MGCD{

RS_THREAD_ATTR pn currentprime;

class Primes:public Crt{
    private:
        static int pr_is_prime_bruteforce(pn n){
            int i,sqrtn;
            sqrtn=(pn)(sqrt((double)n));
            for(i=3;i<=sqrtn;++i)
                if(!(n%i))
                    return 0;
            return 1;
        }

        // vzGG, p. 507; returns 0 if n is composite
        static int pr_fermat(pn n){
            p_set_prime(n);
            return(p_pow(2+((pn)cgal_rs_random())%(n-4),n-1)==1);
        }

        // Solovay-Strassen
        static int pr_ss(pn n){
            pn a,x;
            p_set_prime(n);
            a=1+(pn)cgal_rs_random()%(n-2);
            x=p_div(a,n);
            return(!x||p_pow(a,n>>1)!=x);
        }

        // Miller-Rabin
        static int pr_mr(pn n){
            pn s,d,a,x,r;
            s=1;
            d=(n-1)>>1;
            while(!(d&1)){
                ++s;
                d=d>>1;
            }
            p_set_prime(n);
            a=2+(pn)cgal_rs_random()%(n-4);
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
        static int pr_mrj(pn n){
            pn s,d,a[3],r;//,x;
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

        static pn pr_actual(){
            return currentprime;
        }

    protected:
        static int pr_init(){
            currentprime=PR_MIN;
            return 0;
        }

        static pn pr_next(){
            do{
                currentprime+=2;
            }while(!pr_is_prime(currentprime));
            return currentprime;
        }

}; // class Primes

} // namespace RS_MGCD
} // namespace CGAL

#endif  // CGAL_RS__PRIMES_H

// vim: tabstop=4: softtabstop=4: smarttab: shiftwidth=4: expandtab
