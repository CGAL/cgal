// Copyright (c) 2007 Inria Lorraine (France). All rights reserved.
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

#ifndef CGAL_RS__UGCD_H
#define CGAL_RS__UGCD_H

#include <gmp.h>
#include "primes.h"

// let's assume that 300 is enough for degree 500 gcds
#define CGALRS_MOD_QTY 300

namespace CGAL{
namespace RS_MGCD{

class Ugcd:public Primes{
    public:
        static int ugcd (mpz_t *gcd,mpz_t *Anp,int degA,mpz_t *Bnp,int degB){
            mpz_t *A,*B;
            mpz_t lcgcd,cA,cB;
            mpz_ptr m,bound;
            // dG is initialized to zero only to avoid compiler complaints
            int dA,dB,dG=0,maxd,i,maxA,maxB;
            size_t modsize,modalloc;
            std::vector<CGALRS_PN* > p;
            CGALRS_PN *mA,*mB,*mG,*mod;
            CGALRS_PN lc=0,scaleG;

            if(degB>degA){
                if(!degA){
                    mpz_set_ui(gcd[0],1);
                    return 0;
                }else
                    return ugcd(gcd,Bnp,degB,Anp,degA);
            }
            if(!degB){
                mpz_set_ui(gcd[0],1);
                return 0;
            }

            // initialize the memory
            meminit();

            A=(mpz_t*)malloc((1+degA)*sizeof(mpz_t));
            B=(mpz_t*)malloc((1+degB)*sizeof(mpz_t));
            mpz_init_set(cA,Anp[degA]);
            for(i=0;i<degA;++i)
                mpz_gcd(cA,cA,Anp[i]);
            mpz_init_set(cB,Bnp[degB]);
            for(i=0;i<degB;++i)
                mpz_gcd(cB,cB,Bnp[i]);

            for(i=0;i<=degA;++i){
                mpz_init(A[i]);
                mpz_divexact(A[i],Anp[i],cA);
            }
            for(i=0;i<=degB;++i){
                mpz_init(B[i]);
                mpz_divexact(B[i],Bnp[i],cB);
            }

            // calculate the gcd of the principal coefficients
            mpz_init(lcgcd);
            mpz_gcd(lcgcd,A[degA],B[degB]);

            // find the limit of modular image computation
            maxA=degA;
            for(i=0;i<degA;++i)
                if(mpz_cmpabs(A[i],A[maxA])>0)
                    maxA=i;
            maxB=degB;
            for(i=0;i<degB;++i)
                if(mpz_cmpabs(B[i],B[maxB])>0)
                    maxB=i;
            mpz_pow_ui(cA,A[maxA],2);
            mpz_mul_ui(cA,cA,degA+1);
            mpz_mul_2exp(cA,cA,2*degA);
            mpz_sqrt(cA,cA);
            mpz_pow_ui(cB,B[maxB],2);
            mpz_mul_ui(cB,cB,degB+1);
            mpz_mul_2exp(cB,cB,2*degB);
            mpz_sqrt(cB,cB);
            if(mpz_cmp(cA,cB)<0){
                bound=(mpz_ptr)cA;
                m=(mpz_ptr)cB;
            }else{
                bound=(mpz_ptr)cB;
                m=(mpz_ptr)cA;
            }
            mpz_mul(bound,bound,lcgcd);
            mpz_mul_2exp(bound,bound,1);
            mpz_setbit(bound,0);

            mA=(CGALRS_PN*)palloc((1+degA)*sizeof(CGALRS_PN));
            mB=(CGALRS_PN*)palloc((1+degB)*sizeof(CGALRS_PN));
            maxd=degA;      // we know that degA>=degB
            mG=(CGALRS_PN*)palloc((1+maxd)*sizeof(CGALRS_PN));
            pr_init();
            mpz_set_ui(m,1);
            mod=(CGALRS_PN*)palloc(CGALRS_MOD_QTY*sizeof(CGALRS_PN));
            modalloc=CGALRS_MOD_QTY;
            modsize=0;

            while(mpz_cmp(m,bound)<=0){
                do{
                    p_set_prime(pr_next());
                    dA=pp_from_poly(mA,A,degA);
                    if(dA!=-1){
                        dB=pp_from_poly(mB,B,degB);
                        if(dB!=-1)
                            lc=mpz_fdiv_ui(lcgcd,p_prime());
                        // lc is the image of the principal coefficient
                    }
                }while(dA==-1||dB==-1||!lc
                        ||mpz_divisible_ui_p(A[degA],p_prime())
                        ||mpz_divisible_ui_p(B[degB],p_prime()));
                // now we calculate the gcd mod p_prime
                dG=pp_gcd(mG,mA,degA,mB,degB);
                scaleG=CGALRS_P_DIV(lc,mG[dG]);
                mG[dG]=lc;
                for(i=0;i<dG;++i)
                    mG[i]=p_mul(mG[i],scaleG);
                if(!dG){        // done, we know that the gcd is constant
                    mpz_set_ui(gcd[0],1);
                    dG=0;
                    goto cleanandexit;
                }
                if(dG<maxd){
                    CGALRS_mpz_set_pn(m,p_prime());
                    maxd=dG;
                    p.clear();
                    p.push_back(mG);
                    mod[0]=p_prime();
                    modsize=1;
                    mG=(CGALRS_PN*)palloc((1+maxd)*sizeof(CGALRS_PN));
                    // TODO: clean the  CGALRS_PN* that are in p
                }else{
                    if(dG==maxd){
                        CGALRS_mpz_mul_pn(m,m,p_prime());
                        if(modsize==modalloc){
                            modalloc*=2;
                            mod=(CGALRS_PN*)
                                prealloc(mod,modalloc*sizeof(CGALRS_PN));
                        }
                        mod[modsize]=p_prime();
                        ++modsize;
                        p.push_back(mG);
                        mG=(CGALRS_PN*)palloc((1+maxd)*sizeof(CGALRS_PN));
                    }
                }
            }

            pcra(gcd,mod,p,maxd,modsize);

            cleanandexit:

            CGALRS_PFREE(mA);
            CGALRS_PFREE(mB);
            CGALRS_PFREE(mG);
            CGALRS_PFREE(mod);
            // TODO: clean the CGALRS_PN* that are in p
            for(i=0;i<=degA;++i)
                mpz_clear(A[i]);
            for(i=0;i<=degB;++i)
                mpz_clear(B[i]);
            mpz_clear(m);
            mpz_clear(bound);
            mpz_clear(lcgcd);
            free(A);
            free(B);

            memrelease();

            return dG;
        };

}; // class Ugcd

} // namespace RS_MGCD
} // namespace CGAL

#endif  // CGAL_RS__UGCD_H
